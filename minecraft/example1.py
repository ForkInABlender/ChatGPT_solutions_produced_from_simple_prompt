import cppyy

# Load JVM
cppyy.add_include_path('<path to jvm header files>')
cppyy.add_library_path('<path to jvm dynamic libraries>')
cppyy.load_library('jvm')

# Start JVM
jvm_options = ["-Djava.class.path=<path to minecraft jar file>"]
cppyy.gbl.PyInit_JavaVM(len(jvm_options), jvm_options)

# Import required Java classes
Minecraft = cppyy.gbl.net.minecraft.client.Minecraft
Item = cppyy.gbl.net.minecraft.item.Item
ItemStack = cppyy.gbl.net.minecraft.item.ItemStack
Block = cppyy.gbl.net.minecraft.block.Block
BlockPos = cppyy.gbl.net.minecraft.util.math.BlockPos
EnumFacing = cppyy.gbl.net.minecraft.util.EnumFacing
GuiScreen = cppyy.gbl.net.minecraft.client.gui.GuiScreen

def create_chest():
    mc = Minecraft.getMinecraft()
    inventory = mc.player.inventory.mainInventory
    empty_slots = 0
    for slot in inventory:
        if slot is None:
            empty_slots += 1
    if empty_slots <= 1:
        return  # Player has only one or zero empty slots, no need to create chest
    for i, slot in enumerate(inventory):
        if slot is None:
            # Create chest at player's location
            pos = mc.player.getPosition()
            chest = Block.getBlockFromName("chest")
            mc.world.setBlockState(BlockPos(pos.getX(), pos.getY() + 1, pos.getZ()), chest.getDefaultState())

            # Transfer all items from player's inventory to chest
            for j, item_stack in enumerate(inventory):
                if item_stack is None:
                    continue
                item = item_stack.getItem()
                if item != Item.getItemFromBlock(chest):
                    mc.playerController.windowClick(0, j, 0, 1, mc.player)
                    mc.playerController.windowClick(0, i, 0, 1, mc.player)
                    if empty_slots == 2:
                        break

            # Close chest GUI
            mc.player.closeScreen()
            break

# Test create_chest function
create_chest()
