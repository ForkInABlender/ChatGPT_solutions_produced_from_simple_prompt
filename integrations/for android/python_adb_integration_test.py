# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)

"""

pip install pure-python-adb adb-shell



What this does is allow one to interact with a android device using ADB. This allows one to use python to interact with software running on an android device by its serial id.

This also lets you make use of anything adb normally touches on.

The overrides applied to the shell and more.

I won't go into more details on how as that override can apply even to golang functions integrated also into python. 

https://g.co/bard/share/260ae988f242
"""

from ppadb.client import Client as AdbClient
# Default is "127.0.0.1" and 5037
client = AdbClient(host="127.0.0.1", port=5037)
print(client.version())
print(client.device('serialID-go-here').shell("echo hello world !"))
