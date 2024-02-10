# Dylan Kenneth Eliot & GPT-4-Plugins (Beta Edition)


"""


Now amongst other things, one may load and use minecraft this way. 

On top of extending and overloading classes, you can even wrap them and even overload them.

This is the most affective way to load minecraft into inside of python I have found thus far.

All the main function needs is an access token and a version number. The rest is up to the mod maker from this point onward BEFORE you start the main class.

"""


import jpype
import jpype.imports
import sys

sys_argv=sys.argv[1:]

jvm_args = []
program_args = []
classpath = []

for arg in sys_argv:
    if arg.startswith('-D') or arg.startswith('-X'):
        jvm_args.append(arg)
    elif arg == '-cp':
        classpath.append(sys_argv[sys_argv.index(arg) + 1])
    else:
        program_args.append(arg)

# Convert classpath list to a string
classpath_str = ":".join(classpath)

# Start the JVM
if not jpype.isJVMStarted():
    jpype.startJVM(jpype.getDefaultJVMPath(), *jvm_args, classpath=classpath_str)



import net

net.minecraft.client.main.Main().main(jvm_args+program_args+classpath)


jpype.shutdownJVM()
