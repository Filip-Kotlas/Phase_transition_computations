"""
Converts configs from the old layout (the physical parameters in the top-level sections "problem",
"force_term" and "anisotropy") to the current one, where they are in the section of the model,
model.<type>, in place:

    python3 python/migrate_config.py config/config.json results/OldRun/info/setup/config.json ...

The old force term "zirconium" gives

    "model": {"type": "Zr-Nb_Choud2013",
              "Zr-Nb_Choud2013": {"alpha", "a", "b", "ksi", "A", "m", "theta_0", "d", "T"}}

and "constant" or "radial" gives the model allen_cahn, with "force_term" and, for "constant",
"force_size" instead of "d" and "T". The other sections (solver, init_cond, bound_cond) stay as
they are. Missing old values get the defaults the old solver used. beta was never used by the
solver and is dropped, as are d and T for the allen_cahn model, which does not use them. Files
already in the current layout are left alone.
"""
import json
import sys
from pathlib import Path

# The defaults of the old Parameters::load for the moved values.
OLD_DEFAULTS = {"alpha": 1.0, "a": 1.0, "b": 0.1, "ksi": 0.01, "d": 5e15, "T": 1200.0,
                "A": 0.3, "m": 2, "theta_0": 0.0}
# The sections the model takes over.
OLD_MODEL_SECTIONS = ("problem", "force_term", "anisotropy")


def migrate(config: dict) -> dict:
    """The config in the current layout; the other sections keep their order and content."""
    problem = config.get("problem", {})
    force_term = config.get("force_term", {})
    anisotropy = config.get("anisotropy", {})

    def old(section: dict, key: str):
        return section.get(key, OLD_DEFAULTS[key])

    model_section = {key: old(problem, key) for key in ("alpha", "a", "b", "ksi")}
    model_section.update({key: old(anisotropy, key) for key in ("A", "m", "theta_0")})

    force_type = force_term.get("type", "zirconium")
    if force_type == "zirconium":
        model_type = "Zr-Nb_Choud2013"
        model_section["d"] = old(problem, "d")
        model_section["T"] = old(problem, "T")
    elif force_type in ("constant", "radial"):
        model_type = "allen_cahn"
        model_section["force_term"] = force_type
        if force_type == "constant":
            model_section["force_size"] = force_term.get("size", 0.0)
    else:
        raise ValueError(f"Unknown force term type {force_type!r}.")
    model = {"type": model_type, model_type: model_section}

    # "model" takes the place of the first of the sections it replaces.
    result = {}
    for key, value in config.items():
        if key in OLD_MODEL_SECTIONS:
            result.setdefault("model", model)
        else:
            result[key] = value
    result.setdefault("model", model)
    return result


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for name in sys.argv[1:]:
        path = Path(name)
        config = json.loads(path.read_text(encoding="utf-8"))
        if "model" in config and not any(key in config for key in OLD_MODEL_SECTIONS):
            print(f"{path}: already in the current layout, left alone.")
            continue
        path.write_text(json.dumps(migrate(config), indent=4) + "\n", encoding="utf-8", newline="\n")
        print(f"{path}: converted.")


if __name__ == "__main__":
    main()
