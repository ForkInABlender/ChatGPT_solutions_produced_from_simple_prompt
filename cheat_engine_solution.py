import ctypes

# Load the Cheat Engine DLL
ce = ctypes.WinDLL('cheatengine-x86_64.dll')

# Attach to a process by ID
process_id = 1234
ce.OpenProcess(process_id)

# Get the search value from the user
search_value = int(input("Enter search value: "))

# Search for the value in the process's memory
search_result = ctypes.c_uint64(0)
ce.findWhatWrites(search_value, ctypes.byref(search_result))

# Get the address of the first result
result_address = search_result.value

# Store all the results in a list
results = []
while result_address:
    results.append(result_address)
    result_address = ce.getNextAddress(result_address)

# Modify the values in the process's memory
new_value = int(input("Enter new value: "))
for address in results:
    ce.writeProcessMemory(address, ctypes.byref(new_value), 4, None)
