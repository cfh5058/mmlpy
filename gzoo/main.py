#!/usr/bin/python 
import os,glob,pprint
import numpy as np


def test(clf,datadir=None):
    """Test model"""
    from sklearn.feature_extraction import DictVectorizer
    # Directories
    datadir='/scratch/langmm/gzoo'
    imgtemp=os.path.join(datadir,'images_training_rev1/{}.jpg')
    outfile=os.path.join(datadir,'training_solutions_rev1.csv')
    # Read solutions
    Y,Ykeys=loadout(outfile,idx=slice(50))
    # Get image data
    X,Xkeys = loadimg(imgtemp,Y[:,Ykeys.index('GalaxyID')])
    # Make prediction
    Yp=clf.predict(X)
    # Get error
    rms=np.sqrt(np.mean((Y-Yp)**2))
    return rms

def train(datadir=None):
    """Train model"""
    from sklearn.feature_extraction import DictVectorizer
    from sklearn import svm,feature_selection
    from sklearn import linear_model,ensemble,neighbors
    # Directories
    datadir='/scratch/langmm/gzoo'
    imgtemp=os.path.join(datadir,'images_training_rev1/{}.jpg')
    outfile=os.path.join(datadir,'training_solutions_rev1.csv')
    # Read targets
    Y,Ykeys = loadout(outfile,idx=slice(50,100))
    # Get image data
    X,Xkeys = loadimg(imgtemp,Y[:,Ykeys.index('GalaxyID')])
    # Select features for each target
    # featVal={} ; featKey=np.array(Xkeys)
    # for idx,k in enumerate(Ykeys):
    #     if k=='GalaxyID': continue
    #     y=Y[:,idx]
    #     clf0=feature_selection.SelectKBest(feature_selection.f_regression,k=1)
    #     clf0.fit(X,y)
    #     featVal[k]=clf0.get_support()
    #     print '{}: {}'.format(k,featKey[featVal[k]])
    # Fit
    clf = neighbors.KNeighborsRegressor()
    # clf = ensemble.ExtraTreesRegressor(n_estimators=10, max_features=len(Xkeys),random_state=0)
    # clf = linear_model.LinearRegression()
    # clf = svm.SVR()
    clf.fit(X,Y)
    return clf

# def loadall(meth,datadir=None,**exkw):
#     """Load all jpeg files"""
#     datadir='/scratch/langmm/gzoo'
#     imagdir=os.path.join(datadir,'images_{}_rev1'.format(meth))
#     outfile=os.path.join(datadir,'{}_solutions_rev1.csv'.format(meth))
#     # Loop over files adding color components
#     img_r=[] ; img_g=[] ; img_b=[]
#     flist=glob.glob(imgtemp.format('*'))
#     for f in flist: 
#         ir,ig,ib=loadjpg(f,**exkw)
#         img_r.append(ir)
#         img_g.append(ig)
#         img_b.append(ib)
#     # Flatten filters
#     img_r=np.vstack(img_r)
#     img_g=np.vstack(img_g)
#     img_b=np.vstack(img_b)


def feat_gauss(img):
    """Get gaussian fit parameters"""
    import scipy.optimize as opt
    keys=['b','a','x0','y0','xsig','ysig','rota']
    delkeys=['b','x0','y0','rota']
    Nimg=int(np.sqrt(len(img)))
    X,Y=np.indices((Nimg,Nimg),dtype=float)
    X=X.ravel()
    Y=Y.ravel()
    cw=float(Nimg)/10.
    #p0=[b,a,center_x,center_y,width_x,width_y,rota]
    p0=[min(img),max(img),cw,cw,cw,cw,0.]
    popt, pcov = opt.curve_fit(gauss2D,(X,Y),img,p0=p0,ftol=1.0e5)
    pdct = {k:v for k,v in zip(keys,popt)}
    a=max(pdct['xsig'],pdct['ysig'])
    b=min(pdct['xsig'],pdct['ysig'])
    pdct['ellip']=1.-b/a
    pdct['grit']=np.sqrt(np.mean((img-gauss2D((X,Y),*popt))**2.))
    for k in delkeys: del pdct[k]
    return pdct

def feat_sersic(img):
    """Get sersic fit parameters"""
    import scipy.optimize as opt
    keys=['n','b','a','x0','y0','xsig','ysig','rota']
    Nimg=int(np.sqrt(len(img)))
    X,Y=np.indices((Nimg,Nimg),dtype=float)
    X=X.ravel()
    Y=Y.ravel()
    cw=float(Nimg)/2.
    #p0=[b,a,center_x,center_y,width_x,width_y,rota]
    p0=[2.,min(img),max(img),cw,cw,cw,cw,0.]
    popt, pcov = opt.curve_fit(sersic2D,(X,Y),img,p0=p0)
    pdct = {k:v for k,v in zip(keys,popt)}
    a=max(pdct['xsig'],pdct['ysig'])
    b=min(pdct['xsig'],pdct['ysig'])
    pdct['ellip']=1.-b/a
    pdct['grit']=np.sqrt(np.mean((img-sersic2D((X,Y),*popt))**2.))
    return pdct

def feat_fourier(img):
    """Get fourier tranform parameters"""
    # Nimg=int(np.sqrt(len(img)))
    # X,Y=np.indices((Nimg,Nimg),dtype=float)
    # X=X.ravel()-Nimg/2.
    # Y=Y.ravel()-Nimg/2.
    # R=np.sqrt(X*X+Y*Y)
    # Phi=np.arcsin(Y/R)
    out=np.fft.fft2(img)
    return out

def loadout(fname,idx=None):
    """Routine to load results file"""
    import csv
    from sklearn.feature_extraction import DictVectorizer
    r=csv.DictReader(open(fname,'r'),delimiter=',')
    ydata=[]
    for row in r: ydata.append(np.array([row[k] for k in r.fieldnames],dtype=float))
    if idx: ydata=ydata[idx]
    Y=np.vstack(ydata,)
    names=r.fieldnames
    #print Y.shape,len(names),type(Y),Y.mean()
    return Y,names

def writeout(fname,r,fieldnames=None):
    """Routine to write results file"""
    import csv
    if not fieldnames: fieldnames=getattr(r,'fieldnames',sorted(r[0].keys()))
    f=open(fname,'w')
    w=csv.DictWriter(f,delimiter=',',fieldnames=fieldnames)
    w.writerow(dict((fn,fn) for fn in fieldnames))
    for row in r: w.writerow(row)
    f.close()
    return

def plotjpg(img,fname):
    """Plot image file"""
    import matplotlib.pyplot as plt
    N=int(np.sqrt(len(img[0])))
    for i in range(3):
        plt.subplot(1,3,i+1)
        plt.imshow(img[i].reshape(N,N))
    plt.savefig(fname)
    print '    '+fname
    raise Exception(fname)

def loadimg(ftemp,idList=None,features='gauss'):
    """Load and append all jpeg images"""
    import glob
    from sklearn.feature_extraction import DictVectorizer
    # Get list of files
    if idList is None: flist=sorted(glob.glob(ftemp.format('*')))
    else             : flist=[ftemp.format(int(g)) for g in idList]
    # Define feature generating function
    if   features=='pixels': ffunc=lambda x: x
    elif features=='gauss' : ffunc=lambda x: feat_gauss(x)
    elif features=='sersic': ffunc=lambda x: feat_sersic(x)
    else: raise Exception('Invalid value for features: {}'.format(features))
    # Loop over adding image data
    xdata=[]
    for fimg in flist:
        ir,ig,ib=loadjpg(fimg)
        im=(ir+ig+ib)/3.
        xdata.append(ffunc(im))
    # Vectorize
    if features=='pixels': 
        X=np.vstack(xdata)
        names=np.arange(len(xdata))
    else:
        xv=DictVectorizer()
        X=xv.fit_transform(xdata).toarray()
        names=xv.get_feature_names()
    # Return data
    return X,names

def loadjpg(fname,method='PIL',size=(50,50),plot=False):
    """Routine to load jpeg file"""
    if method=='PIL':
        from PIL import Image
        import Image
        im = Image.open(fname).resize(size)
        px = np.array(im.getdata()).reshape((size[1],size[0],3))
        px = [px[:,:,i].ravel() for i in range(3)]
        if plot: plotjpg(px,fname.replace('.jpg','_rgb.png'))
        return px


def sersic2D(xy,n=2.,b=0.,a=1.,center_x=0.,center_y=0.,width_x=1.,width_y=1.,rota=0.):
    """Returns amplitude of 2D sersic profile at specified location"""
    if n==0: return np.zeros(len(xy[0]))
    x,y=xy
    rcen_x = center_x*np.cos(rota) - center_y*np.sin(rota)
    rcen_y = center_x*np.sin(rota) + center_y*np.cos(rota)
    xp = x*np.cos(rota) - y*np.sin(rota)
    yp = x*np.sin(rota) + y*np.cos(rota)
    g = b + a*np.exp(-(((rcen_x-xp)/width_x)**(1./float(n))+
                       ((rcen_y-yp)/width_y)**(1./float(n)))/2.)
    return g.ravel()

def gauss2D(xy,b=0.,a=1.,center_x=0.,center_y=0.,width_x=1.,width_y=1.,rota=0.):
    """Returns amplitude of 2D gaussian at specified location"""
    x,y=xy
    rcen_x = center_x*np.cos(rota) - center_y*np.sin(rota)
    rcen_y = center_x*np.sin(rota) + center_y*np.cos(rota)
    xp = x*np.cos(rota) - y*np.sin(rota)
    yp = x*np.sin(rota) + y*np.cos(rota)
    g = b + a*np.exp(-(((rcen_x-xp)/width_x)**2+
                       ((rcen_y-yp)/width_y)**2)/2.)
    return g.ravel()

