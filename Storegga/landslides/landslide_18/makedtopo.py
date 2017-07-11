"""
Module to combine fixed_grids output into one dtopo file

by Jihwan Kim

"""
def dtopo_tt1():
   from numpy import mod,zeros

   fid=open('_output/fixed_grids.data','r')
   for i in range(8):
      line=fid.readline()

   line=fid.readline().split()
   ng=int(line[2])

   fid.close()

   infile="fort.fg01_xxxx"
   dir="_output/"

   odir='/media/jihwan/EXT_NGI/Computational_Results/Bing_FVM/remolding/slope_1d/dtopo'
   outfile=odir+"/slide_2k.tt1"

   infile1=dir+"fort.fg01_0001"
   fid= open(infile1,'r')
   fout=open(outfile,'w')

   line=fid.readline().split()
   time=float(line[0])

   line=fid.readline().split()
   mx=int(line[0])

   line=fid.readline().split()
   my=int(line[0])

   line=fid.readline().split()
   xlow=float(line[0])

   line=fid.readline().split()
   ylow=float(line[0])

   line=fid.readline().split()
   xhi=float(line[0])

   line=fid.readline().split()
   yhi=float(line[0])

   line=fid.readline()
   line=fid.readline()

   h_orig=zeros(mx*my)
   k=0
   dx=(xhi-xlow)/(mx-1)
   dy=(yhi-ylow)/(my-1)
   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h_orig[k]=float(line[0])
         xx=xlow+i*dx
         yy=yhi-j*dy
         k+=1
         fout.write("%20.10e %20.10e %20.10e %20.10e \n" % (time, xx, yy, 0.e0))

   fid.close()

   infile1=infile
   dz=zeros(mx*my)

   for k in range(2,ng+1):
      for ipos in range(14,10,-1):
         idigit=mod(k,10)
         infile1=infile1[0:ipos-1]+str(idigit)+infile1[ipos:14]
         k=k/10

      fid1= open(dir+infile1,'r')
      print fid1
      line=fid1.readline().split()
      tt=float(line[0])

      for i in range(8):
         line=fid1.readline().split()
  
      k=0  
      for j in range(my):
         for i in range(mx):
            line1=fid1.readline().split()
            hnew=float(line1[0])
            dz[i+(my-j-1)*mx]=hnew-h_orig[k]
            k+=1
      
      for j in range(my):
         for i in range(mx):    
            xx=xlow+i*dx
            yy=yhi-j*dy
            k=i+j*mx
      
            fout.write("%20.10e %20.10e %20.10e %20.10e \n" % (tt, xx, yy, dz[k]))

      fid1.close()

   fout.close()

def dtopo_tt3(indir,outfile):
   """
   Generated tt3 file from fixedgrid output
   """
   from numpy import mod,zeros

   fid_data = indir+'/fixed_grids.data'

   fid=open(fid_data,'r')
   for i in range(8):
      line=fid.readline()

   line=fid.readline().split()
   t0  =float(line[0])
   tf  =float(line[1])
   mt  =int(line[2])
   xlow=float(line[3])
   xup =float(line[4])
   ylow=float(line[5])
   yup =float(line[6])
   mx  =int(line[7])
   my  =int(line[8])

   fid.close()

   dx =(xup-xlow)/(mx-1)
   dy =(yup-ylow)/(my-1)
   dt =(tf-t0)/(mt-1)

   dir = indir

   infile="fort.fg01_xxxx"

   infile1=dir+"/fort.fg01_0001"
   fid= open(infile1,'r')
   fout=open(outfile,'w')

   fout.write("%10d \n" % mx)
   fout.write("%10d \n" % my)
   fout.write("%10d \n" % (mt-1))
   fout.write("%20.10e \n" % xlow)
   fout.write("%20.10e \n" % ylow)
   fout.write("%20.10e \n" % (t0+dt))
   fout.write("%20.10e \n" % dx)
   fout.write("%20.10e \n" % dy)
   fout.write("%20.10e \n" % dt)

   for i in range(9):
      line=fid.readline()

   h_orig=zeros(mx*my)
   k=0
   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h_orig[k]=float(line[0])
         k+=1

   fid.close()

   infile1=infile
   dz=zeros(mx*my)

   for k in range(2,mt+1):
      for ipos in range(14,10,-1):
         idigit=mod(k,10)
         infile1=infile1[0:ipos-1]+str(idigit)+infile1[ipos:14]
         k=k/10

      fid1= open(dir+'/'+infile1,'r')
      print fid1

      for i in range(9):
         line=fid1.readline().split()
  
      k=0  
      for j in range(my):
         for i in range(mx):
            line1=fid1.readline().split()
            hnew=float(line1[0])
            dz[i+(my-j-1)*mx]=hnew-h_orig[k]
            k+=1
      
      for j in range(my):
         for i in range(mx):    
            k=i+j*mx
      
            fout.write("%20.10e " % dz[k])

         fout.write("\n")

      fid1.close()

   fout.close()

def dtopo_tt3_v2(indir,outfile):
   """
   Generated tt3 file from fort.q output
   """
   from numpy import mod,zeros

   fid_data = indir+'/claw.data'

   fid=open(fid_data,'r')
   for i in range(15):
      line=fid.readline()

   line=fid.readline().split()
   t0  =float(line[0])

   for i in range(2):
      line=fid.readline()

   line=fid.readline().split()
   mt  =int(line[0])

   line=fid.readline().split()
   tf  =float(line[0])

   fid.close()

   dt =(tf-t0)/mt

   dir = indir

   infile="fort.qxxxx"

   infile1=dir+"/fort.q0000"
   fid= open(infile1,'r')
   fout=open(outfile,'w')

   for i in range(2):
      line=fid.readline()

   line=fid.readline().split()
   mx = int(line[0])
   line=fid.readline().split()
   my = int(line[0])
   line=fid.readline().split()
   xlow = float(line[0])
   line=fid.readline().split()
   ylow = float(line[0])
   line=fid.readline().split()
   dx = float(line[0])
   line=fid.readline().split()
   dy = float(line[0])

   line=fid.readline()

   fout.write("%10d \n" % mx)
   fout.write("%10d \n" % my)
   fout.write("%10d \n" % mt)
   fout.write("%20.10e \n" % xlow)
   fout.write("%20.10e \n" % ylow)
   fout.write("%20.10e \n" % (t0+dt))
   fout.write("%20.10e \n" % dx)
   fout.write("%20.10e \n" % dy)
   fout.write("%20.10e \n" % dt)

   h_orig=zeros(mx*my)
   k=0
   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h_orig[k]=float(line[0])
         k+=1
      line=fid.readline()

   fid.close()

   infile1=infile
   dz=zeros(mx*my)

   for k in range(1,mt+1):
      for ipos in range(10,6,-1):
         idigit=mod(k,10)
         infile1=infile1[0:ipos-1]+str(idigit)+infile1[ipos:10]
         k=k/10

      fid1= open(dir+'/'+infile1,'r')
      print fid1

      for i in range(9):
         line=fid1.readline().split()
  
      k=0  
      for j in range(my):
         for i in range(mx):
            line1=fid1.readline().split()
            hnew=float(line1[0])
            dz[i+(my-j-1)*mx]=hnew-h_orig[k]
            k+=1
         line=fid1.readline()
      
      for j in range(my):
         for i in range(mx):    
            k=i+j*mx
      
            fout.write("%20.10e " % dz[k])

         fout.write("\n")

      fid1.close()

   fout.close()

def slide_eta_tt3(outfile):
   """
   Generated tt3 file from fort.q0000 output
   """
   from numpy import mod,zeros

   sea_level = 0.

   infile="fort.q0000"
   dir="_output/"
   infile1=dir+infile

   fid= open(infile1,'r')
   fout=open(outfile,'w')

   grid=fid.readline()
   amr=fid.readline()
   mx=int(fid.readline().split()[0])
   my=int(fid.readline().split()[0])
   xlow=float(fid.readline().split()[0])
   ylow=float(fid.readline().split()[0])
   dx=float(fid.readline().split()[0])
   dy=fid.readline()
   nodata=fid.readline()
   nodata=-9999

   fout.write("%10d \n" % mx)
   fout.write("%10d \n" % my)
   fout.write("%20.10e \n" % xlow)
   fout.write("%20.10e \n" % ylow)
   fout.write("%20.10e \n" % dx)
   fout.write("%20.10e \n" % nodata)

   h=zeros((mx,my))
   eta=zeros((mx,my))

   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h[i,j]=float(line[0])
         eta[i,j]=float(line[6]) + sea_level
      fid.readline()
   fid.close()


   for j in range(my):
      for i in range(mx):
         fout.write("%20.10e " % eta[i,my-j-1])

      fout.write("\n")

   fout.close()

def slide_h_tt1(outfile):
   """
   Generated tt3 file from fort.q0000 output
   """
   from numpy import mod,zeros

   infile="fort.q0000"
   dir="_output/"
   infile1=dir+infile

   fid= open(infile1,'r')
   fout=open(outfile,'w')

   grid=fid.readline()
   amr=fid.readline()
   mx=int(fid.readline().split()[0])
   my=int(fid.readline().split()[0])
   xlow=float(fid.readline().split()[0])
   ylow=float(fid.readline().split()[0])
   dx=float(fid.readline().split()[0])
   dy=float(fid.readline().split()[0])
   nodata=fid.readline()
   nodata=-99999

   h_orig=zeros((mx,my))

   for j in range(my):
      for i in range(mx):
         line=fid.readline().split()
         h_orig[i,j]=float(line[0])
      fid.readline()
   fid.close()

   for j in range(my):
      y = ylow+(my-j-.5)*dy
      for i in range(mx):
          x = xlow+(i+.5)*dx
          fout.write("%20.10e %20.10e %20.10e \n" % (x,y, h_orig[i,my-j-1]))
   fout.close()

if __name__=='__main__':
    outfile="dtopo_storegga_10kto500.tt3"
    #dtopo_tt3(indir,outfile)
    outfile="slide_eta_0m.tt3"
    slide_eta_tt3(outfile)
    #outfile="slide_h.tt1"
    #slide_h_tt1(outfile)
    indir = '_output'
    outfile="dtopo_storegga_10kto500.tt3"
    #dtopo_tt3_v2(indir,outfile)

