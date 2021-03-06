import numpy as np
import ROOT as rt
#import matplotlib.pyplot as plt
#from matplotlib.patches import RegularPolygon, Polygon
import csv
import argparse as arg

parser = arg.ArgumentParser(description='Train Binary BDT')
parser.add_argument('-l', '--layer', dest='layer', type=str, default='1', help="layer number")
parser.add_argument('-v', '--version', dest='version', type=str, default='v16', help="layer number")
parser.add_argument('-z', '--zside', dest='side', type=str, default='p', help="layer number")

args = parser.parse_args()
layer = args.layer
version = args.version
zside = args.side

if version == 'v15':
    rootfile = 'geantoutput_v15_28.root'
    csvfile = 'wafer_v15.csv'
else:
    rootfile = 'geantoutput_v16_2_11.root'
    csvfile = 'wafer_v16.csv'

if len(layer) == 1:
    layer_name = '0' + layer
else:
    layer_name = layer

wafer_size = 166.44/10
wafer_saperation = 1/10

r = 0.5 * (wafer_saperation + wafer_size)
R = 2 * r / np.sqrt(3)
dy = 0.75 * R

u = []
v = []
wafer_thick = [] 
wafer_type = []
wafer_orient = []


with open(csvfile, mode = 'r') as file:
    csvfile = csv.reader(file)
    for lines in csvfile:
        if lines[1] != layer: continue
        u.append(float(lines[2]))
        v.append(float(lines[3]))
        wafer_thick.append(lines[5])
        wafer_type.append(lines[6])
        wafer_orient.append(float(lines[7]))

s3 = np.sqrt(3)
partial = {}
partial['full'] = np.array([[s3/2, 1/2], [0, 1], [-s3/2, 1/2], [-s3/2, -1/2], [0, -1], [s3/2, -1/2], [s3/2, 1/2]])
partial['a'] = np.array([[0, 1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [0, 1]])
partial['d'] = np.array([[-s3/4, 3/4], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [s3/4, -3/4], [-s3/4, 3/4]])
partial['g'] = np.array([[s3/4, 3/4], [0,1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [s3/4, -3/4], [s3/4, 3/4]])
partial['b'] = np.array([[s3/2, 1/2], [0,1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [s3/2, 1/2]])
partial['c'] = np.array([[-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [-s3/2, 1/2]])
partial['gm'] = np.array([[1/(2*s3), 5/6], [0,1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [1/(2*s3), -5/6], [1/(2*s3), 5/6]])
partial['dm'] = np.array([[-1/s3, 2/3], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [1/(2*s3), -5/6], [-1/s3, 2/3]])

#u = [0, 0, 0, 1, 1, 1, -1, -1, -1]
#v = [0, 1, -1, 0, 1, -1, 0, 1, -1]
print(np.shape(partial['a']))

iu = np.array(u)
iv = np.array(v)

x = -1 * ((-2 * iu + iv) * r)
y = 2 * iv * dy

tfile = rt.TFile(rootfile, 'READ')
if zside == 'n':
    hXYF = tfile.Get('prodEE/gXYhitsF0_layer_' + layer_name)
    hXYC1 = tfile.Get('prodEE/gXYhitsCN0_layer_' + layer_name)
    hXYC2 = tfile.Get('prodEE/gXYhitsCK0_layer_' + layer_name)
    #hXYB = tfile.Get('prodEE/hXYhitsB_layer_' + layer_name)
else:
    hXYF = tfile.Get('prodEE/gXYhitsF1_layer_' + layer_name)
    hXYC1 = tfile.Get('prodEE/gXYhitsCN1_layer_' + layer_name)
    hXYC2 = tfile.Get('prodEE/gXYhitsCK1_layer_' + layer_name) 

print(hXYF)
#print(hXYF.Integral())
if zside == 'n':
    hXYF.GetXaxis().SetTitle('x (cm)')
    hXYF.SetTitle('Hits in XY for layer ' + layer + ' (-z side)')
elif zside == 'p':
    hXYF.GetXaxis().SetTitle('-x (cm)')
    hXYF.SetTitle('Hits in XY for layer ' + layer + ' (+z side)')
elif zside == 'pp':
    hXYF.GetXaxis().SetTitle('x (cm)')
    hXYF.SetTitle('Hits in XY for layer ' + layer + ' (+z side)')
    for i in range(hXYF.GetN()): hXYF.SetPointX(i, -1*hXYF.GetPointX(i))
    for i in range(hXYC1.GetN()): hXYC1.SetPointX(i, -1*hXYC1.GetPointX(i))
    for i in range(hXYC2.GetN()): hXYC2.SetPointX(i, -1*hXYC2.GetPointX(i)) 
    #hXYF.Scale(-1, 'x')
    #hXYC1.Scale(-1, 'x')
    #hXYC2.Scale(-1, 'x')
hXYF.GetYaxis().SetTitle('y (cm)')
hXYF.GetYaxis().SetTitleOffset(1.4)
hXYF.GetXaxis().SetLimits(-200.0, 200.0)
hXYF.GetYaxis().SetRangeUser(-200.0, 200.0)

hXYF.SetLineColor(rt.kRed)
hXYF.SetMarkerColor(rt.kRed)
hXYF.SetFillColor(rt.kRed)
hXYF.SetMarkerStyle(8)
hXYF.SetMarkerSize(0.2)

hXYC1.SetLineColor(rt.kGreen + 1)
hXYC1.SetMarkerColor(rt.kGreen + 1)
hXYC1.SetFillColor(rt.kGreen + 1)
hXYC1.SetMarkerStyle(8)
hXYC1.SetMarkerSize(0.2)

hXYC2.SetLineColor(rt.kMagenta)
hXYC2.SetMarkerColor(rt.kMagenta)
hXYC2.SetFillColor(rt.kMagenta)
hXYC2.SetMarkerStyle(8)
hXYC2.SetMarkerSize(0.2)

leg = rt.TLegend(0.78,0.15,0.99,0.3)
leg.AddEntry(hXYF,"Si width 120 #mum","lf")
leg.AddEntry(hXYC1,"Si width 200 #mum","lf")
leg.AddEntry(hXYC2,"Si width 300 #mum","lf")
#leg.AddEntry(hXYhitsB[il],"CEH Sci","lf")



#print(x,y)
c1 = rt.TCanvas('c1', 'Layer ' + layer, 800, 800, 800, 800)
c1.cd()
graphs = []

hXYF.Draw('ap')
hXYC1.Draw('p,SAME')
hXYC2.Draw('p,SAME')
leg.Draw()

# Add some coloured hexagons
for i in range(len(x)):
    thetad = 60*wafer_orient[i]
    theta = np.radians(thetad)
    #theta = 0
    rot_mat = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    edges = R * np.matmul(partial[wafer_type[i]], rot_mat)
    shift_edges = edges + np.array([[x[i], y[i]]])
    if zside == 'pp':
        ix = -1 * np.array(shift_edges[:, 0])
    else:
        ix = np.array(shift_edges[:, 0])
    iy = np.array(shift_edges[:, 1])
    n = np.shape(ix)[0]
    #print(n)
    graphs.append(rt.TGraph(n, ix, iy))
    graphs[-1].Draw('l, SAME')
    #print(np.transpose(np.matmul(rot_mat, np.transpose(partial['c']))))

# Also add scatter points in hexagon centres
c1.Update()
c1.Show()
c1.Print('over_l' + layer + '_' + version + '_2_' + zside + '.png')
#print(u, v, wafer_orient, wafer_type, wafer_thick)

input()
#print(r, R, dy)
