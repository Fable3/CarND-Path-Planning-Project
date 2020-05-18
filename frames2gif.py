import imageio
import os
import cv2

image_folder="frames_sel/"
image_list = sorted([os.path.join(image_folder, image_file) for image_file in os.listdir(image_folder)])
                        
with imageio.get_writer('passing_sample.gif', mode='I', duration=0.5) as writer:
    for filename in image_list:
        image = imageio.imread(filename)
        image=cv2.resize(image, (640,360))
        writer.append_data(image)
