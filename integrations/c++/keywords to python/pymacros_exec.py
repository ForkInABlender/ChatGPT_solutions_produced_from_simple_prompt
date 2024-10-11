# Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )

"""
This enabled defining macros, a basepath, and code block that can be executed.

Similar to the c++ preprocessor, minus a few macros like `#progma`, but with the exception of conditional handling.

`#ifdef` checks if a macro is defined conditionally. `#ifndef` does similar for undefined's.
This was done to handle conditional checks like `#if !defined{macro} || {optional condition}` for the purposes of mimicking how the c++ preprocessors do it before your code is normally compiled down.

This allows for templating the includes, if's, else's,endif's, ifdef's, ifndef's and `#if !defined`'s declarative control of expression.

This, similarly enables preprocessor conditional definitions of code for ctypes' struct's & union's.

"""

def preprocess_code(code, macros, base_path='.'):
    lines = code.strip().split('\n')
    output_code = []
    conditional_block = False

    for line in lines:
        line_strip = line.strip()
        if line_strip.startswith('#ifdef'):
            # Enhanced to handle expressions after macro definition
            parts = line_strip.split(maxsplit=1)
            if len(parts) > 1:
                try:
                    # Evaluate expression if additional conditions exist beyond the macro name
                    conditional_block = parts[1] not in macros or not eval(parts[1], {}, {'True': True, 'False': False, **macros})
                except Exception as e:
                    print(f"Error evaluating '#ifdef' condition '{parts[1]}': {e}")
                    conditional_block = True  # Default to true if evaluation fails
            else:
                conditional_block = parts[1] not in macros
        elif line_strip.startswith('#ifndef'):
            # Enhanced to handle expressions after macro definition
            parts = line_strip.split(maxsplit=1)
            if len(parts) > 1:
                try:
                    # Evaluate expression if additional conditions exist beyond the macro name
                    conditional_block = parts[1] in macros or eval(parts[1], {}, {'True': True, 'False': False, **macros})
                except Exception as e:
                    conditional_block = False  # Default to false if evaluation fails
                    pass
            else:
                conditional_block = parts[1] in macros
        elif line_strip.startswith('#if'):
            expression = line_strip[3:].strip()
            try:
                conditional_block = not eval(expression, {}, {'True': True, 'False': False, **macros})
            except Exception as e:
                print(f"Error evaluating expression '{expression}': {e}")
                conditional_block = True  # Default to true if evaluation fails
        elif line_strip.startswith('#else'):
            conditional_block = not conditional_block
        elif line_strip.startswith('#endif'):
            conditional_block = False
        elif line_strip.startswith('#define'):
            parts = line_strip.split(maxsplit=2)
            if len(parts) == 3:
                macros[parts[1]] = eval(parts[2], {}, macros) if any(c.isdigit() for c in parts[2]) else parts[2]
        elif line_strip.startswith('#undef'):
            parts = line_strip.split()
            if len(parts) == 2:
                macros.pop(parts[1], None)
        elif line_strip.startswith('#include'):
            parts = line_strip.split(maxsplit=1)
            if len(parts) == 2:
                filename = parts[1].strip('"').strip("<>").strip()
                file_path = base_path+filename
                try:
                    with open(file_path, 'r') as file:
                        included_code = file.read()
                    output_code.append(preprocess_code(included_code, macros, base_path=base_path))
                except FileNotFoundError:
                    raise Exception(f"File not found: {file_path}")
        elif not conditional_block:
            for key, value in macros.items():
                line = line.replace(str(key), str(value))
            output_code.append(line)

    return '\n'.join(output_code)
