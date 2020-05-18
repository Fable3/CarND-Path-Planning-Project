import cv2
print(cv2.__version__)
vidcap = cv2.VideoCapture('passing_sample.mkv')
success,image = vidcap.read()
count = 0
while success:
  cv2.imwrite("frames/frame%04d.png" % count, image)
  success,image = vidcap.read()
  count += 1
  
import imageio
import os
image_folder="frames/"
image_list = sorted([os.path.join(image_folder, image_file) for image_file in os.listdir(image_folder)])
                        
with imageio.get_writer('passing_sample.gif', mode='I') as writer:
    for filename in image_list:
        image = imageio.imread(filename)
        writer.append_data(image)
        