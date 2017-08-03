#
#
#

# import my stuff
import os, shutil, sys, getopt
from subprocess import call

# parse input arguments
if len(sys.argv) < 6:
    print 'image_generate.py <inputimage> <new_width> <xdim> <outputname> <n_mc>'
    sys.exit()

inputimage = sys.argv[1];
width = sys.argv[2];
xdim = sys.argv[3];
outputname = sys.argv[4];
n_mc = sys.argv[5];

Sigma = [0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256]

print "------------------------------------------" 
print " name of file  : %s" % (inputimage)
print " new width     : %s" % (width)
print " xdim          : %s" % (xdim)
print " output name   : %s" % (outputname)
print " n_mc          : %s" % (n_mc)
print "------------------------------------------" 

outputdirectory = 'data/%s' % (outputname)

# create folder
if not os.path.exists(outputdirectory):
    os.makedirs(outputdirectory)

# generate solution
filename = "solution.bin";
os.system("./util_dlib_image_to_vec --filename_in=%s --new_width=%s --xdim=%s --type=0 --noise=%s --filename_out=%s/%s" % (inputimage, width, xdim, 0.0,outputdirectory,filename) )
#os.system("./util_dlib_image_to_vec --filename_in=%s --xdim=%s --type=0 --noise=%s --filename_out=%s/%s" % (inputimage, xdim, 0.0,outputdirectory,filename) )
os.remove("%s/%s.info" % (outputdirectory, filename));

for n in range(int(n_mc)):
    for idSigma in range(len(Sigma)):
        filename = "%s_id%s_idSigma%s.bin" % (outputname,n,idSigma)
        print "generating %s/%s" % (outputdirectory,filename)
        os.system("./util_dlib_image_to_vec --filename_in=%s --new_width=%s --xdim=%s --type=0 --noise=%s --filename_out=%s/%s" % (inputimage, width, xdim, Sigma[idSigma],outputdirectory,filename) )
        os.remove("%s/%s.info" % (outputdirectory, filename));


