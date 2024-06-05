To use the `inject_so.py` script and target the current bash shell process, run the following command:

```
python3.10 inject_so.py $(cat /proc/self/status | grep '^Pid:' | awk '{print $2}') /path/to/libexample.so
```

This command targets the bash shell process from which it is executed. The script will run in the background, and it will inject the shared object into the bash shell process as 
 long as the shell remains active.

https://github.com/ancat/gremlin is the github to reach into for the inject_so.py script. You may have to rewrite it for python3.10 syntax.
It is recommended that such run as a background process for affective single shell hacks. 
