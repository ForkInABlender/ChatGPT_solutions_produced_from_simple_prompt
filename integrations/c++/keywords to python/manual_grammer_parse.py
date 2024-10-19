# Dylan Kenneth Eliot & GPT-4o ( Canvas; Alpha Edition )

"""
This parses reinterpret_cast<type>(expr) as seen in c++ code.

Yes, this uses regular expression to enforce proper ctypes utilization.

From here, it is a matter of mapping keyword usage over syntax layout and templating so it makes use of each keyword correctly.

Because this is part of the project to remove the need for a compiler, this may also be useful for moving c++ projects from c++ to pure python code instead.


"""

import ctypes
import re

def reinterpret_cast(original_value, target_type, original_type=None):
    """
    Reinterpret the bits of an original value as a different ctypes type.
    Equivalent to reinterpret_cast<target_type>(original_value) in C++.

    :param original_value: The value to reinterpret.
    :param target_type: The ctypes type to cast to (e.g., ctypes.c_float).
    :param original_type: The original ctypes type (e.g., ctypes.c_int). If None, infer from original_value.
    :return: The reinterpreted value.
    """
    # Infer original type if not provided
    if original_type is None:
        if isinstance(original_value, ctypes._SimpleCData):
            original_type = type(original_value)
        else:
            raise ValueError("original_type must be specified for non-ctypes values")

    # Ensure the original value is an instance of the original type
    if not isinstance(original_value, original_type):
        original_value = original_type(original_value)

    # Cast the original pointer to the target type pointer
    target_ptr = ctypes.cast(ctypes.pointer(original_value), ctypes.POINTER(target_type))
    return target_ptr.contents.value

def reinterpret_cast_from_string(input_string):
    """
    Parse and return reinterpret_cast operations from a given string.

    :param input_string: The input string containing reinterpret_cast operations.
    :return: A string containing the modified code that can be manually executed.
    """
    # Replace all occurrences of reinterpret_cast with a callable Python function
    def replace_reinterpret_cast(match):
        target_type_str = match.group(1)
        original_value_str = match.group(2).strip()
        return f"reinterpret_cast({original_value_str}, ctypes.{target_type_str})"

    # Update the pattern to match reinterpret_cast expressions
    pattern = r"reinterpret_cast<([a-zA-Z_][a-zA-Z0-9_]*)>\((.+?)\)"
    modified_code = re.sub(pattern, replace_reinterpret_cast, input_string)

    # Correcting indentation errors by stripping leading/trailing whitespace from each line
    corrected_code_lines = [line.strip() for line in modified_code.splitlines() if line.strip()]
    corrected_code = "\n".join(corrected_code_lines)

    # Return the modified code to be executed manually
    return corrected_code
