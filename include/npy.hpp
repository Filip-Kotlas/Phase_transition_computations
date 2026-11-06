#pragma once

// Minimal reader and writer of the NumPy .npy format (version 1.0), C-ordered arrays of float or double.
// The files are little-endian, so the code assumes a little-endian machine (x86, ARM in its usual mode).

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace npy {

// Type descriptor of T in the header, e.g. "<f8" for double.
template< typename T >
std::string descriptor()
{
    static_assert( std::is_same< T, float >::value || std::is_same< T, double >::value,
                   "Only float and double are supported." );
    return std::is_same< T, double >::value ? "<f8" : "<f4";
}

inline const unsigned char* magic()
{
    static const unsigned char bytes[ 6 ] = { 0x93, 'N', 'U', 'M', 'P', 'Y' };
    return bytes;
}

// Writes the header of a C-ordered array of the given shape. The raw data follow it directly.
inline void write_header( std::ostream& out, const std::string& descr, const std::vector< std::size_t >& shape )
{
    std::string shape_text = "(";
    for( std::size_t i = 0; i < shape.size(); i++ )
        shape_text += std::to_string( shape[ i ] ) + ( shape.size() == 1 ? "," : ( i + 1 < shape.size() ? ", " : "" ) );
    shape_text += ")";

    std::string header = "{'descr': '" + descr + "', 'fortran_order': False, 'shape': " + shape_text + ", }";
    // 6 bytes of magic string, 2 of version, 2 of header length and the header ending with '\n'
    // have to take a multiple of 64 bytes together.
    const std::size_t total = 10 + header.size() + 1;
    header.append( ( 64 - total % 64 ) % 64, ' ' );
    header += '\n';

    const std::uint16_t length = static_cast< std::uint16_t >( header.size() );
    const unsigned char version_and_length[ 4 ] = { 1, 0,
                                                    static_cast< unsigned char >( length & 0xFF ),
                                                    static_cast< unsigned char >( length >> 8 ) };
    out.write( reinterpret_cast< const char* >( magic() ), 6 );
    out.write( reinterpret_cast< const char* >( version_and_length ), 4 );
    out.write( header.data(), static_cast< std::streamsize >( header.size() ) );
}

// Reads the header written by NumPy (or write_header). The stream is left at the first byte of the data.
inline void read_header( std::istream& in, std::string& descr, bool& fortran_order, std::vector< std::size_t >& shape )
{
    unsigned char start[ 8 ];
    in.read( reinterpret_cast< char* >( start ), 8 );
    if( ! in || std::memcmp( start, magic(), 6 ) != 0 )
        throw std::runtime_error( "The file is not in the .npy format." );

    std::size_t header_length = 0;
    if( start[ 6 ] == 1 ) {
        unsigned char length[ 2 ];
        in.read( reinterpret_cast< char* >( length ), 2 );
        header_length = length[ 0 ] | ( static_cast< std::size_t >( length[ 1 ] ) << 8 );
    }
    else if( start[ 6 ] == 2 || start[ 6 ] == 3 ) {
        unsigned char length[ 4 ];
        in.read( reinterpret_cast< char* >( length ), 4 );
        for( int i = 0; i < 4; i++ )
            header_length |= static_cast< std::size_t >( length[ i ] ) << ( 8 * i );
    }
    else
        throw std::runtime_error( "Unsupported version of the .npy format." );

    std::string header( header_length, ' ' );
    in.read( &header[ 0 ], static_cast< std::streamsize >( header_length ) );
    if( ! in )
        throw std::runtime_error( "The .npy header is truncated." );

    // The header is a Python dictionary: {'descr': '<f8', 'fortran_order': False, 'shape': (2, 3, 4), }
    auto value_start = [ & ]( const std::string& key ) {
        std::size_t pos = header.find( key );
        if( pos == std::string::npos )
            throw std::runtime_error( "The .npy header has no " + key + "." );
        pos = header.find( ':', pos );
        if( pos == std::string::npos )
            throw std::runtime_error( "The .npy header is malformed." );
        return pos + 1;
    };

    std::size_t pos = value_start( "'descr'" );
    const std::size_t quote_begin = header.find( '\'', pos );
    const std::size_t quote_end = quote_begin == std::string::npos ? quote_begin : header.find( '\'', quote_begin + 1 );
    if( quote_end == std::string::npos )
        throw std::runtime_error( "The .npy header is malformed." );
    descr = header.substr( quote_begin + 1, quote_end - quote_begin - 1 );

    pos = header.find_first_not_of( ' ', value_start( "'fortran_order'" ) );
    fortran_order = pos != std::string::npos && header.compare( pos, 4, "True" ) == 0;

    pos = value_start( "'shape'" );
    const std::size_t open = header.find( '(', pos );
    const std::size_t close = header.find( ')', open );
    if( open == std::string::npos || close == std::string::npos )
        throw std::runtime_error( "The .npy header is malformed." );
    shape.clear();
    bool in_number = false;
    for( std::size_t i = open + 1; i < close; i++ ) {
        if( header[ i ] >= '0' && header[ i ] <= '9' ) {
            if( ! in_number )
                shape.push_back( 0 );
            shape.back() = shape.back() * 10 + static_cast< std::size_t >( header[ i ] - '0' );
            in_number = true;
        }
        else
            in_number = false;
    }
}

}  // namespace npy
