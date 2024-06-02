import re


def parse_argument(arg_str: str):
    arg_type, arg_name = re.match(r"([\w\s]+\s+\*{0,})\s{0,}(\w+)$", arg_str).groups()
    return arg_type.strip(), arg_name


def parse_header_file(filename):
    funcs = []
    enums = dict()
    func_str = ""
    enum_str = ""

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if enum_str or line.startswith("typedef enum"):
                enum_str += line

                if enum_str.endswith(";"):
                    re_match = re.match(r"typedef enum (\w+) \{(.+)\} (\w+);", enum_str)
                    assert re_match.group(1) == re_match.group(3)
                    enum_vals = [x.strip() for x in re_match.group(2).split(",")]

                    enums[re_match.group(3)] = enum_vals
                    enum_str = ""

            # only support void return type
            elif func_str or line.startswith("void"):
                func_str += line

                if func_str.endswith(";"):
                    re_match = re.match(r"void ([\w_]+)\((.+)\);", func_str)
                    func_args = [parse_argument(x) for x in re_match.group(2).split(",")]

                    funcs.append((re_match.group(1), func_args))
                    func_str = ""

    return funcs, enums


def generate_py_func(func_name: str, func_args: list[tuple[str, str]], enums: dict[str, list[str]]):
    arg_pattern_lookup = {
        "const char *": "y#",
        "char *": "y*",
        "int": "i",
        "float": "f",
        "double": "d",
    }

    c_vars = []
    c_func_args = []
    arg_parse_pattern = ""

    for arg_type, arg_name in func_args:
        if arg_type in enums:
            arg_type = "int"

        # read-only buffer
        if arg_type.startswith("const ") and arg_type.endswith("*"):
            arg_type = "const char *"
            c_vars.append((arg_type, arg_name))
            c_vars.append(("Py_ssize_t", f"{arg_name}_size"))
            c_func_args.append(arg_name)

        elif arg_type.endswith("*"):
            arg_type = "char *"
            c_vars.append(("Py_buffer", arg_name))
            c_func_args.append(f"{arg_name}.buf")

        elif arg_type in ("int", "float", "double"):
            c_vars.append((arg_type, arg_name))
            c_func_args.append(arg_name)

        else:
            raise ValueError(f"{func_name} - Unsupported {arg_type=} for {arg_name=}")

        arg_parse_pattern += arg_pattern_lookup[arg_type]

    lines = [f"static PyObject *py_{func_name}(PyObject *self, PyObject *args) {{"]

    # declare local variables
    for c_type, c_name in c_vars:
        lines.append(f"  {c_type} {c_name};")
    lines.append("")

    # parse input arguments
    lines.append(
        f'  if (!PyArg_ParseTuple(args, "{arg_parse_pattern}", {", ".join("&" + var_name for _, var_name in c_vars)}))'
    )
    lines.append("    return NULL;")
    lines.append("")
    lines.append("  Py_BEGIN_ALLOW_THREADS;")
    lines.append(f'  {func_name}({", ".join(c_func_args)});')
    lines.append("  Py_END_ALLOW_THREADS;")
    lines.append("")
    for var_type, var_name in c_vars:
        if var_type == "Py_buffer":
            lines.append(f"  PyBuffer_Release(&{var_name});")
    lines.append("")
    lines.append("  Py_INCREF(Py_None);")
    lines.append("  return Py_None;")
    lines.append("}")
    return lines


def generate_extension_source(header_file):
    funcs, enums = parse_header_file(header_file)

    lines = [
        '#include "img_proc.h"',
        "#define PY_SSIZE_T_CLEAN",
        "#include <Python.h>",
        "",
    ]

    for func_name, func_args in funcs:
        lines.extend(generate_py_func(func_name, func_args, enums))
        lines.append("")

    lines.append("static PyMethodDef ImgProcMethods[] = {")
    for func_name, _ in funcs:
        lines.append(f'  {{"{func_name}", py_{func_name}, METH_VARARGS, ""}},')
    lines.append("  {NULL, NULL, 0, NULL},")  # sentinel
    lines.append("};")
    lines.append("")
    lines.append("static struct PyModuleDef img_proc_module = {")
    lines.append('    PyModuleDef_HEAD_INIT, "img_proc", "module docstring", -1, ImgProcMethods,')
    lines.append("};")
    lines.append("")
    lines.append("PyMODINIT_FUNC PyInit_img_proc(void) { return PyModule_Create(&img_proc_module); }")
    return lines


if __name__ == "__main__":
    from pathlib import Path

    CURRENT_DIR = Path(__file__).parent
    filename = Path(__file__).name

    lines = generate_extension_source(CURRENT_DIR / "img_proc.h")
    with open(CURRENT_DIR / "py_img_proc.c", "w") as f:
        f.write(f"// This file is auto-generated from {filename}\n")
        f.write("\n".join(lines))
        f.write("\n")
