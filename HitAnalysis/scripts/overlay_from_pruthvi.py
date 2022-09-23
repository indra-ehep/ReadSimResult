import numpy as np
import ROOT as rt
#import matplotlib.pyplot as plt
#from matplotlib.patches import RegularPolygon, Polygon
import csv
import argparse as arg

def GetCartesian(r, theta):
  x = r*rt.TMath.Cos(theta*rt.TMath.Pi()/180.)
  y = r*rt.TMath.Sin(theta*rt.TMath.Pi()/180.)
  return x,y;

def GetPolar(x, y):
  r = rt.TMath.Sqrt(x*x + y*y)
  theta =  180.0*rt.TMath.ATan2(y, x)/rt.TMath.Pi()
  if theta<0.0:
    theta = (180. - rt.TMath.Abs(theta))  + 180.    
  return r,theta;


parser = arg.ArgumentParser(description='Train Binary BDT')
parser.add_argument('-l', '--layer', dest='layer', type=str, default='28', help="layer number")
parser.add_argument('-v', '--version', dest='version', type=str, default='v17', help="version number")
parser.add_argument('-z', '--zside', dest='side', type=str, default='n', help="zside number")

args = parser.parse_args()
layer = args.layer
version = args.version
zside = args.side

if version == 'v15':
    rootfile = 'geantoutput_v15_28.root'
    csvfile = 'wafer_v15.csv'
elif version == 'v16':
    rootfile = 'geantoutput_v16_2_11.root'
    csvfile = 'wafer_v16.csv'
elif version == 'v17':
    rootfile = 'geantoutput_D92.root'
    csvfile = 'wafer_v17.csv'



wafer_size = 166.44/10
wafer_saperation = 1/10

r = 0.5 * (wafer_saperation + wafer_size)
R = 2 * r / np.sqrt(3)
dy = 0.75 * R

cassette_shifts = {}
layer_types = {}

"""
with open(csvfile, mode = 'r') as file:
    csvfile = csv.reader(file)
    for lines in csvfile:
        if lines[1] != layer: continue
        u.append(float(lines[2]))
        v.append(float(lines[3]))
        wafer_thick.append(lines[5])
        wafer_type.append(lines[6])
        wafer_orient.append(float(lines[7]))
"""
thick_map = {'h120': 'Fine', 'l200': 'Coarse1', 'l300': 'Coarse2'}
type_thick = {'h120': 'HD' , 'l200': 'LD', 'l300': 'LD'}
flatfile = open('flatfile.txt', 'r')
data = flatfile.read().split('\n')
for line in data[:47]:
    lines = line.split()
    layer_types[lines[0]] = lines[1]
    cassette_shifts[lines[0]] = [float(shift) for shift in lines[2:]]

sci_sensor_file = open('Validation/HGCalValidation/data/scintillatorV0.txt','r')
sci_data = sci_sensor_file.read().split('\n')
#print(sci_data)
for line in sci_data[:14]:
    lines = line.split()
    # if len(lines)==10 :
    #     print(lines)
    
s3 = np.sqrt(3)
points = [[0, -1], [s3/2, -1/2], [s3/2, 1/2], [0, 1], [-s3/2, 1/2], [-s3/2, -1/2], [-s3/4, -3/4], [s3/4, 3/4], [-0.6855, -0.6043], [-0.1806, -0.8937], [0.6855, -0.6043], [0.6855, 0.6043], [0.1806, 0.8957], [-s3/2, 0.2914]]
partial = {}
partial['full'] = np.array([[s3/2, 1/2], [0, 1], [-s3/2, 1/2], [-s3/2, -1/2], [0, -1], [s3/2, -1/2], [s3/2, 1/2]])
partial['a'] = np.array([[0, 1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [0, 1]])
partial['d'] = np.array([[-s3/4, 3/4], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [s3/4, -3/4], [-s3/4, 3/4]])
partial['g'] = np.array([[s3/4, 3/4], [0,1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [s3/4, -3/4], [s3/4, 3/4]])
partial['b'] = np.array([[s3/2, 1/2], [0,1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [s3/2, 1/2]])
partial['c'] = np.array([[-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [-s3/2, 1/2]])
partial['gm'] = np.array([[1/(2*s3), 5/6], [0,1], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [1/(2*s3), -5/6], [1/(2*s3), 5/6]])
partial['dm'] = np.array([[-1/s3, 2/3], [-s3/2, 1/2], [-s3/2, -1/2], [0,-1], [1/(2*s3), -5/6], [-1/s3, 2/3]])
partial['LD1'] = np.array([points[4], points[5], points[0], points[1], points[4]])
partial['LD2'] = np.array([points[1], points[2], points[3], points[4], points[1]])
partial['LD3'] = np.array([points[0], points[1], points[2], points[7], points[6], points[0]])
partial['LD4'] = np.array([points[6], points[7], points[3], points[4], points[5], points[6]])
partial['LD5'] = np.array([points[0], points[1], points[2], points[3], points[5], points[0]])
partial['LD6'] = np.array([points[5], points[3], points[4], points[5]])
partial['HD1'] = np.array([points[0], points[10], points[13], points[5], points[0]])
partial['HD2'] = np.array([points[10], points[1], points[2], points[3], points[4], points[13], points[10]])
partial['HD3'] = np.array([points[0], points[1], points[2], points[11], points[9], points[0]])
partial['HD4'] = np.array([points[8], points[12], points[3], points[4], points[5], points[8]])
partial['HD5'] = np.array([points[0], points[1], points[2], points[12], points[8], points[0]])

tfile = rt.TFile(rootfile, 'READ')

for lay in range(47,48):
    layer = str(lay)
    if len(layer) == 1:
        layer_name = '0' + layer
    else:
        layer_name = layer

    
    layer_type = layer_types[layer]
    cassette_shift = cassette_shifts[layer]
    
    for ishift in range(0,len(cassette_shift)) :
      cassette_shift[ishift] = float(cassette_shift[ishift])/10.0
    print(cassette_shift, " cm")
    
    #print(layer_type)
    
    u = []
    v = []
    wafer_thick = [] 
    wafer_type = []
    wafer_orient = []
    cassette = []
    
    for line in data[47:13589]:
        lines = line.split()
        #print('kill', lines)
        #print(line)
        if lines[0] != layer: continue
        u.append(float(lines[6]))
        v.append(float(lines[7]))
        wafer_thick.append(thick_map[lines[2]])
        if lines[1] == '0':
            wafer_type.append('full')
        else:
            wafer_type.append(type_thick[lines[2]]+ lines[1])
        wafer_orient.append(float(lines[5]))
        cassette.append(int(lines[8]))

    #print(u, v, wafer_thick,  wafer_type,  wafer_orient, cassette)



    #u = [0, 0, 0, 1, 1, 1, -1, -1, -1]
    #v = [0, 1, -1, 0, 1, -1, 0, 1, -1]
    

    iu = np.array(u)
    iv = np.array(v)

    if layer_type == '1':
        x = 1 * ((-2 * iu + iv) * r)
    else:
        x = -1 * ((-2 * iu + iv) * r)
    y = 2 * iv * dy

    if layer_type == '3':
        y = y - R

    if layer_type == '2':
        y = y + R

    print(lay, len(wafer_orient), len(u), len(x), len(y), len(cassette))
    if int(layer) <=26:
        prod = 'prodEE'
        prod1 = 'prodEE'
    else:
        prod = 'prodHEF'
        prod1 = 'prodHEB'
        
    if zside == 'n':
        hXYF = tfile.Get(prod + '/grXYhitsF0_layer_' + layer_name)
        hXYC1 = tfile.Get(prod + '/grXYhitsCN0_layer_' + layer_name)
        hXYC2 = tfile.Get(prod + '/grXYhitsCK0_layer_' + layer_name)
        hXYB1 = tfile.Get(prod1 + '/grXYhitsB0_layer_' + layer_name)
        if version == 'v17':
            hXYAR = tfile.Get(prod + '/grXYhitsAR0_layer_' + layer_name)
            #hXYB = tfile.Get('prodEE/hXYhitsB_layer_' + layer_name)
    else:
        hXYF = tfile.Get(prod + '/grXYhitsF1_layer_' + layer_name)
        hXYC1 = tfile.Get(prod + '/grXYhitsCN1_layer_' + layer_name)
        hXYC2 = tfile.Get(prod + '/grXYhitsCK1_layer_' + layer_name)
        hXYB1 = tfile.Get(prod1 + '/grXYhitsB1_layer_' + layer_name)
        if version == 'v17':
            hXYAR = tfile.Get(prod + '/grXYhitsAR1_layer_' + layer_name)

    print(hXYF)
    #print(hXYF.Integral())
    if zside == 'n':
        hXYC2.GetXaxis().SetTitle('x (cm)')
        hXYC2.SetTitle('Hits in XY for layer ' + layer + ' (-z side)')
    elif zside == 'p':
        hXYC2.GetXaxis().SetTitle('-x (cm)')
        hXYC2.SetTitle('Hits in XY for layer ' + layer + ' (+z side)')
        for i in range(hXYF.GetN()): hXYF.SetPointX(i, -1*hXYF.GetPointX(i))
        for i in range(hXYC1.GetN()): hXYC1.SetPointX(i, -1*hXYC1.GetPointX(i))
        for i in range(hXYC2.GetN()): hXYC2.SetPointX(i, -1*hXYC2.GetPointX(i))
        for i in range(hXYB1.GetN()): hXYB1.SetPointX(i, -1*hXYB1.GetPointX(i))
        for i in range(hXYAR.GetN()): hXYAR.SetPointX(i, -1*hXYAR.GetPointX(i)) 
    elif zside == 'pp':
        hXYF.GetXaxis().SetTitle('x (cm)')
        hXYF.SetTitle('Hits in XY for layer ' + layer + ' (+z side)') 
        #hXYF.Scale(-1, 'x')
        #hXYC1.Scale(-1, 'x')
        #hXYC2.Scale(-1, 'x')
    hXYF.GetYaxis().SetTitle('y (cm)')
    hXYF.GetYaxis().SetTitleOffset(1.4)
    hXYF.GetXaxis().SetLimits(-300.0, 300.0)
    hXYF.GetYaxis().SetRangeUser(-300.0, 300.0)

    hXYC2.GetYaxis().SetTitle('y (cm)')
    hXYC2.GetYaxis().SetTitleOffset(1.4)
    hXYC2.GetXaxis().SetLimits(-300.0, 300.0)
    hXYC2.GetYaxis().SetRangeUser(-300.0, 300.0)
    
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

    if int(layer) >= 34:
        hXYB1.SetLineColor(rt.kBlue)
        hXYB1.SetMarkerColor(rt.kBlue)
        hXYB1.SetFillColor(rt.kBlue)
        hXYB1.SetMarkerStyle(8)
        hXYB1.SetMarkerSize(0.1)
        hXYB1.GetXaxis().SetLimits(-300.0, 300.0)
        hXYB1.GetYaxis().SetRangeUser(-300.0, 300.0)

    if version == 'v17':
        hXYAR.SetLineColor(rt.kBlack)
        hXYAR.SetMarkerColor(rt.kBlack)
        hXYAR.SetFillColor(rt.kBlack)
        hXYAR.SetMarkerStyle(8)
        hXYAR.SetMarkerSize(0.2)

    leg = rt.TLegend(0.78,0.15,0.99,0.3)
    leg.AddEntry(hXYF,"Si width 120 #mum","lf")
    leg.AddEntry(hXYC1,"Si width 200 #mum","lf")
    leg.AddEntry(hXYC2,"Si width 300 #mum","lf")
    leg.AddEntry(hXYB1,"Scintillator","lf")
    #leg.AddEntry(hXYhitsB[il],"CEH Sci","lf")



    #print(x,y)
    c1 = rt.TCanvas('c1', 'Layer ' + layer, 600, 600, 600, 600)
    c1.cd()
    graphs = []
    text = []
    #text = rt.TPaveText(-200, -200, 200, 200)
    if int(layer) <= 33:
        hXYC2.Draw('ap')
        hXYC1.Draw('p,SAME')
        hXYF.Draw('p,SAME')
    else:
        hXYC2.Draw('ap')
        hXYC1.Draw('p,SAME')
        hXYF.Draw('p,SAME')
        hXYB1.Draw('p,SAME')
    if version == 'v17':
        hXYAR.Draw('p,SAME')
        leg.Draw()

    # Add some coloured hexagons
    for i in range(len(x)):
        thetad = 60*wafer_orient[i]
        theta = np.radians(thetad)
        #theta = 0
        rot_mat = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
        edges = R * np.matmul(partial[wafer_type[i]], rot_mat)
        shift_edges = edges + np.array([[x[i], y[i]]]) + np.array([cassette_shift[2*cassette[i]-2], cassette_shift[2*cassette[i]-1]])
        #print(layer_type == '1')
        if (bool(zside == 'pp') ^ bool(layer_type == '1')):
            ix = -1 * np.array(shift_edges[:, 0])
            if version == 'v17':
                text.append(rt.TPaveText(-1 * x[i] - 3, y[i] - 3, -1 * x[i] + 3, y[i] + 3, 'NB'))
        else:
            ix = np.array(shift_edges[:, 0])
            text.append(rt.TPaveText(x[i] - 3, y[i] - 3, x[i] + 3, y[i] + 3, 'NB'))
        iy = np.array(shift_edges[:, 1])
        n = np.shape(ix)[0]
        #text[-1].AddText('(' + str(u[i]) + ', ' + str(v[i]) + ')')
        graphs.append(rt.TGraph(n, ix, iy))
        graphs[-1].Draw('l, SAME')
        if version == 'v17':
            text[-1].AddText(str(int(wafer_orient[i])))
            text[-1].SetFillColorAlpha(0, 0)
            text[-1].Draw('SAME')
            #print(np.transpose(np.matmul(rot_mat, np.transpose(partial['c']))))
            
    # loop over the lines of scintillator sensor file (aka flatfile of scintillator)
    for line in sci_data[2:437]:
        lines = line.split()
        if len(lines)==10 :
            sc_layer, ring, rmin, rmax = lines[0], int(lines[1]), float(lines[2]),float(lines[3])
            sipmsize, hex1, hex2, hex3 = float(lines[4]), lines[5], lines[6], lines[7]
            hex4, scitype = lines[8], lines[9]
            hexcomb = str(bin(int(hex1,16))[2:].zfill(8)) + str(bin(int(hex2,16))[2:].zfill(8)) + str(bin(int(hex3,16))[2:].zfill(8)) + str(bin(int(hex4,16))[2:].zfill(8))
            # To convert to cm from mm input values
            rmax /= 10.0
            rmin /= 10.0
            iseg = 0
            if int(sc_layer)==lay :
                icassette = 0
                #print(lines, hexcomb, rmin, rmax)
                for isect in range(0,3):
                    for iphi in range(0,len(hexcomb)):
                        #print("isect %i, iphi : %i, icassette %i"%(isect,iphi,icassette))
                        phimin = isect*120. + iphi*1.25
                        phimax = isect*120. + (iphi+1)*1.25
                        nrmin, nrmax = rmin, rmax
                        
                        if iseg%24==0 :
                            icassette += 1

                        # ToDo : For the moment interpreting shift values of flat files in mm for overlay plot of scintillator.
                        xshift, yshift = float(cassette_shift[2*(icassette-1)]), float(cassette_shift[2*(icassette-1)+1])

                        
                        x1, y1 = GetCartesian(rmin, phimin)
                        x2, y2 = GetCartesian(rmin, phimax)
                        x3, y3 = GetCartesian(rmax, phimin)
                        x4, y4 = GetCartesian(rmax, phimax)

                        x1, y1 = x1+xshift, y1+yshift
                        x2, y2 = x2+xshift, y2+yshift
                        x3, y3 = x3+xshift, y3+yshift
                        x4, y4 = x4+xshift, y4+yshift
                                                
                        nrmin, phimin = GetPolar(x1, y1)
                        nrmax, phimax = GetPolar(x4, y4)
                        
                        if int(hexcomb[iphi])==1:
                          #print("Arc1 : ",offsetX, offsetY, nrmax, phimin,phimax)
                          arc1 = rt.TArc(0.0, 0.0, nrmax, phimin,phimax)
                          arc1.SetFillStyle(0)
                          arc1.SetNoEdges()
                          arc1.SetLineWidth(1)
                          graphs.append(arc1)
                          graphs[-1].Draw('l, SAME')
                          
                          l1 = rt.TLine(x1,y1,x3,y3)
                          l1.SetLineWidth(1)
                          graphs.append(l1)
                          graphs[-1].Draw('l, SAME')

                          arc2 = rt.TArc(0.0, 0.0, nrmin, phimin,phimax)
                          arc2.SetFillStyle(0)
                          arc2.SetNoEdges()
                          arc2.SetLineWidth(1)
                          graphs.append(arc2)
                          graphs[-1].Draw('l, SAME')

                          l2 = rt.TLine(x2,y2,x4,y4)
                          l2.SetLineWidth(1)
                          graphs.append(l2)
                          graphs[-1].Draw('l, SAME')
                          
                        iseg += 1;
                        
    # Also add scatter points in hexagon centres
    print('(' + str(u[i]) + ', ' + str(v[i]) + ')\n' + str(wafer_orient[i]))
    c1.Update()
    #c1.Show()
    #if lay == 1:
    #    c1.Print('over_lay' + version + '_' + zside + '.pdf[', 'pdf')
    c1.Print('over_lay_l' + layer_name + '_'  + version + 'old1_' + zside + '.png', 'png')
    c1.Print('over_lay_l' + layer_name + '_'  + version + 'old1_' + zside + '.pdf', 'pdf')
#c1.Print('over_lay' + version + '_' + zside + '.pdf]', 'pdf')
#print(u, v, wafer_orient, wafer_type, wafer_thick)

#input()
#print(r, R, dy)

