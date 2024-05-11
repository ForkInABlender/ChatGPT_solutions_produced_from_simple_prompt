#Dylan Kenneth Eliot & GPT-4 ( Alpha Edition )


"""
The script turns only an image into a frame index, coordinates and color list.
The overall purpose of this script is for usage with pybrain3 and numpy+numba cuda simulation
 for further processing with the custom MAML Ai model.


Because all unit, assert, and mock testing are done by hand manually,
 one can be sure it isn't the automation process causing halt in development.
Later it will be used to generate images using pybrain3, numpy, numba, brian2, and rdkit with smile groups.


"""




from PIL import Image
import numpy as np

def generate_updates_from_image(image_path, i):
    img = Image.open(image_path)
    img_array = np.array(img)
    
    updates = []
    frame_index = i  # Assuming there is only one frame
    
    # Iterate over each pixel in the image
    for x in range(img_array.shape[0]):
        for y in range(img_array.shape[1]):
            color = img_array[x, y].tolist()  # Get the color of the pixel
            update = (frame_index, x, y, color)
            updates.append(update)
    
    return updates

# Example usage
image_path = 'path_to_your_image_representing_frame0.png'
updates = generate_updates_from_image(image_path, 0)
print(updates)
