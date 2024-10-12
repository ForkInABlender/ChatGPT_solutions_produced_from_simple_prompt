# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

This enables c++ comments around, above, in and through a function can be added to a singular help string for python.

This relocates the string for proper documantation keeping and code compatibility.

"""

import re

def cpp_comments_to_docstrings(code):
    # Split the code into lines for easier processing
    lines = code.split('\n')
    transformed_lines = []
    
    # This will hold the comments that will be used in the docstring
    current_docstring = []
    
    # Track if we're in a block comment
    in_block_comment = False
    
    for line in lines:
        stripped_line = line.strip()

        # Handle starting of a block comment
        if stripped_line.startswith('/*'):
            in_block_comment = True
            # Extract the block comment
            comment = stripped_line[2:].strip()
            current_docstring.append(comment)
            if '*/' in stripped_line:  # If the block comment closes on the same line
                in_block_comment = False
                current_docstring[-1] = current_docstring[-1].split('*/')[0].strip()
            continue
            
        # If we are inside a block comment, continue collecting comments
        if in_block_comment:
            if '*/' in stripped_line:  # Closing the block comment
                in_block_comment = False
                current_docstring.append(stripped_line.split('*/')[0].strip())
            else:
                current_docstring.append(stripped_line)
            continue
        
        # Detect function definitions
        func_def_match = re.match(r'def\s+(?:/\*.*?\*/\s*)?(\w+)\s*\((.*?)\):', stripped_line)
        if func_def_match:
            func_name = func_def_match.group(1)
            params = func_def_match.group(2)
            
            # Extract inline comments in parameters
            param_comments = re.findall(r'/\*\s*(.*?)\s*\*/', params)
            cleaned_params = re.sub(r'/\*\s*.*?\s*\*/', '', params).strip()

            # Prepare docstring content
            docstring_lines = []
            if current_docstring:
                docstring_lines.extend(current_docstring)
                current_docstring = []
            if param_comments:
                docstring_lines.extend(param_comments)

            # Start the function definition
            transformed_lines.append(f'def {func_name}({cleaned_params}):')
            transformed_lines.append('    """')
            transformed_lines.extend(f'    {line}' for line in docstring_lines)
            transformed_lines.append('    """')

            # Now handle the function body
            continue
        
        # Clean up C++ comments from the line
        # Remove inline block comments
        cleaned_line = re.sub(r'/\*\s*.*?\s*\*/', '', stripped_line)
        # Remove single-line comments
        cleaned_line = re.sub(r'//.*', '', cleaned_line)

        # Append the cleaned line if it's not empty, maintaining indentation
        if cleaned_line.strip():
            transformed_lines.append(line.replace(stripped_line, cleaned_line))  # Keep the original line with indentation

    return '\n'.join(transformed_lines)
