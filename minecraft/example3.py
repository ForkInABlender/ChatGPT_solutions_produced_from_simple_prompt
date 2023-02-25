import cppyy

# Load the JVM from source code using cppyy
cppyy.include("jni.h")
cppyy.include("jvmti.h")
cppyy.include("jvm.h")
cppyy.include("jni_md.h")
cppyy.add_include_path("/usr/lib/jvm/java-11-openjdk-amd64/include")
cppyy.add_include_path("/usr/lib/jvm/java-11-openjdk-amd64/include/linux")
cppyy.load_library("jvm", "/usr/lib/jvm/java-11-openjdk-amd64/lib/server")

# Get the Minecraft instance from the JVM
jvm = cppyy.gbl.JavaVM()
jvm.AttachCurrentThread()
jclass = jvm.FindClass("net/minecraft/client/Minecraft")
minecraft = jvm.CallStaticObjectMethod(jclass, jvm.GetStaticMethodID(jclass, "getInstance", "()Lnet/minecraft/client/Minecraft;"))

# Get the player's inventory
player = minecraft.player
inventory = player.inventory

# Check if there's only one empty slot in the inventory
empty_slots = 0
for i in range(36):
    if inventory.getStack(i).isEmpty():
        empty_slots += 1
if empty_slots == 1:
    # Create a chest at the player's position
    x, y, z = player.posX, player.posY, player.posZ
    player.world.setBlockState(jvm.JInt(x), jvm.JInt(y), jvm.JInt(z), jvm.getStaticField(jvm.FindClass("net/minecraft/init/Blocks"), "CHEST", "Lnet/minecraft/block/Block;").getDefaultState())

# Detach from the JVM
jvm.DetachCurrentThread()
