include 'imager.g'
include 'viewer.g'

# default arguments
imagetype:="observed"
msname := "F"
mode := "mfs";
weighting := "briggs";
stokes := "I";
select := '';
npix := 256;
cell := '1arcsec';
image_viewer := 'kvis';
fieldid := 1;
spwid := 1;
chanmode := 'channel';
nchan := 1;
chanstart := 1;
chanstep := 1;
img_nchan := 1;
img_chanstart := 1;
img_chanstep := 1;
wprojplanes := 0;
padding := 1.0;
fitsname := '';
filter_bmaj := '';
filter_bmin := '';
filter_bpa := '';
phasecenter := '';

# parse command line
for( a in argv )
{
  print 'arg: ',a;
  if( a =='DATA' )
    imagetype:="observed";
  else if( a =='MODEL_DATA' )
    imagetype:="model";
  else if( a =='CORRECTED_DATA' )
    imagetype:="corrected";
  else if( a =='RESIDUAL' )
    imagetype:="residual";
  else if( a =~ s/^ms=// )
    msname := a;
  else if( a =~ s/^fits=// )
    fitsname := a;
  else if( a =~ s/^mode=// )
    mode := a;
  else if( a =~ s/^weight=// )
    weighting := a;
  else if( a =~ s/^stokes=// )
    stokes := a;
  else if( a =~ s/^npix=// )
    npix := as_integer(a);
  else if( a =~ s/^cellsize=// )
    cell := a;
  else if( a =~ s/^viewer=// )
    image_viewer := a;
  else if( a =~ s/^field=// )
    fieldid := as_integer(a);
  else if( a =~ s/^spwid=// )
    spwid := as_integer(a);
  else if( a =~ s/^wprojplanes=// )
    wprojplanes := as_integer(a);
  else if( a =~ s/^padding=// )
    padding := as_float(a);
  else if( a =~ s/^filter_bmaj=// )
    filter_bmaj := a;
  else if( a =~ s/^filter_bmin=// )
    filter_bmin := a;
  else if( a =~ s/^filter_bpa=// )
    filter_bpa := a;
  else if( a =~ s/^phasecenter=// )
    phasecenter := a;
  else if( a =~ s/^chanmode=// )
    chanmode := a;
  else if( a =~ s/^nchan=// )
    nchan := as_integer(a);
  else if( a =~ s/^chanstart=// )
    chanstart := as_integer(a);
  else if( a =~ s/^chanstep=// )
    chanstep := as_integer(a);
  else if( a =~ s/^img_nchan=// )
    img_nchan := as_integer(a);
  else if( a =~ s/^img_chanstart=// )
    img_chanstart := as_integer(a);
  else if( a =~ s/^img_chanstep=// )
    img_chanstep := as_integer(a);
  else if( a =~ s/^select=// )
    select := a;
}
if( select != '' )
  select := spaste('( ',select,' ) && ANTENNA1 != ANTENNA2');
else
  select := 'ANTENNA1 != ANTENNA2';
print "Selection string: ",select;

# setup the imager
myimager:=imager(msname)
myimager.setdata(mode=chanmode,
             fieldid=fieldid,
             spwid=spwid,
             nchan=nchan,
             start=chanstart,
             step=chanstep,
             msselect=select,
             async=F);

myimager.setimage(nx=npix,ny=npix,cellx=cell,celly=cell,
  phasecenter=phasecenter,
  stokes=stokes,mode=mode,
  fieldid=fieldid,spwid=spwid,
  nchan=img_nchan,start=img_chanstart,step=img_chanstep);

if( weighting != 'default' )
  myimager.weight(weighting); 
  
if( filter_bmaj != '' )
  myimager.filter('gaussian',filter_bmaj,filter_bmin,filter_bpa);
  
if( wprojplanes > 0 )
  myimager.setoptions(ftmachine='wproject', wprojplanes=wprojplanes,cache=500000000,padding=padding);
else
  myimager.setoptions(cache=500000000,padding=padding);


imgname := msname
imgname =~ s/\..*//;
imgname =~ s/.*\///;
imgname := spaste(imgname,".",imagetype,"-",stokes,"-",mode,spaste(img_nchan));
imgfile := spaste(imgname,".img");

# make the image
myimager.makeimage(type=imagetype,image=imgfile,async=F);


myimager.done()

# convert to FITS
im := image(imgfile);
if( fitsname == '' )
  fitsname := spaste(imgname,'.fits');
im.tofits(fitsname,overwrite=T);
im.done();
shell(spaste('rm -fr ',imgfile));
print "\n\n--------- wrote FITS image: ",fitsname," ---------\n";

# run Karma
if( image_viewer != '')
{
  cmd := paste(image_viewer,fitsname);
  shell(cmd);
}
exit
