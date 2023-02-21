import ctypes

# Load the GameConqueror DLL
gc = ctypes.CDLL('libgameconqueror.so')

# Attach to a process by ID
process_id = 1234
gc.gc_attach_process(process_id)

# Get the search value from the user
search_value = int(input("Enter search value: "))

# Search for the value in the process's memory
search_result = ctypes.c_void_p(0)
gc.gc_find_all(search_value, ctypes.byref(search_result))

# Get the address of the first result
result_address = search_result.value

# Store all the results in a list
results = []
while result_address:
    results.append(result_address)
    result_address = gc.gc_get_next_address(result_address)

# Modify the values in the process's memory
new_value = int(input("Enter new value: "))
for address in results:
    gc.gc_edit(address, new_value, 4)
