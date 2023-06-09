import cv2
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

def contour_to_segments(contours):
    x0 = []
    y0 = []
    x1 = []
    y1 = []
    for i in range(len(contours)):
        contour = np.squeeze(contours[i])
        for j in range(len(contour)//2):
            # print(contour)
            try:
                x0.append(contour[2*j][0])
                y0.append(contour[2*j][1])
                x1.append(contour[2*j+1][0])
                y1.append(contour[2*j+1][1])
            except:
                pass
    
    return x0, y0, x1, y1

dir = 'mask/stroma/'
print(os.listdir(dir))
for mask in os.listdir(dir):
    name = mask.split('.')[0]
    print(name)
    img = cv2.imread(dir+mask)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    contours, hierarchy = cv2.findContours(img, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    x0, y0, x1, y1 = contour_to_segments(contours)
    df = pd.DataFrame({'x0':x0, 'y0':y0, 'x1':x1, 'y1':y1})
    df.to_csv('h_contour/stroma_csv/'+name+'.csv', index=False)

