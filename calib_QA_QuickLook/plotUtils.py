from matplotlib import pyplot as plt
import numpy as np

import lsst.afw.cameraGeom.utils as cameraGeomUtils
import lsst.afw.display as afwDisplay

def plotHists(expId, butler, dataType, histColor=None, detectorOrder=None, saturationList=None,
              fileUsedDir=None):
    fig, axes = plt.subplots(3,3, figsize=(13, 13), constrained_layout=True)
    histKwargs = dict(histtype='stepfilled', alpha=0.8, bins=4000, fc=histColor, ec="k")
    camera = butler.get("camera")
    saturationList = [1e12]*9
    for ccd in camera:
        amps = ccd.getAmplifiers()
        for amp in amps:
            saturationList[ccd.getId()] = min(saturationList[ccd.getId()], amp.getSaturation())
    print("saturationList = ", saturationList)
    print("")
    print("det    min        max      median    mean    stddev      var")
    formatStr = "{:2d} {:9.2f} {:10.2f} {:8.2f} {:8.2f} {:8.2f} {:9.2f}"
    dataArr1DList = []
    nonSatArrList = []
    xmin, xmax = 1e12, -1e12
    for iCcd, ax in enumerate(axes.flatten()):
        detNum = detectorOrder[iCcd]
        dataId = {'expId': expId, 'detector': detNum}
        dataImage = butler.get(dataType, dataId=dataId)
        dataArr1D = dataImage.maskedImage.image.array.ravel()
        dataArr1DList.append(dataArr1D)
        nonSatArr = dataArr1D[np.abs(dataArr1D) < saturationList[iCcd]]
        nonSatArrList.append(nonSatArr)
        xmin = min(xmin, np.percentile(nonSatArr, 0.5))
        xmax = max(xmax, np.percentile(nonSatArr, 99.5))
    for iCcd, ax in enumerate(axes.flatten()):
        detNum = detectorOrder[iCcd]
        detName = camera[detNum].getName()
        dataArr1D = dataArr1DList[iCcd]
        nData = len(dataArr1D)
        nanmin, nanmax = np.nanmin(dataArr1D), np.nanmax(dataArr1D)
        datamin, datamax = np.min(dataArr1D), np.max(dataArr1D)
        if (datamin != nanmin or datamax != nanmax):
            print("NOTE: there are non-finite numbers in the data array")
        nanmedian, nanmean = np.nanmedian(dataArr1D), np.nanmean(dataArr1D)
        nanstd, nanvar = np.nanstd(dataArr1D), np.nanvar(dataArr1D)
        if iCcd == 0 and nanstd < 1.0:
            formatStr = "{:2d} {:9.3f} {:10.3f} {:8.4f} {:8.4f} {:8.4f} {:9.4f}"
        print(formatStr.format(iCcd, nanmin, nanmax, nanmedian, nanmean, nanstd, nanvar))
        nonSatArr = nonSatArrList[iCcd]
        ax.hist(nonSatArr, label="Det: {} ({})".format(detNum, detName), **histKwargs)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, 600000)
        ax.set_title("Det: {} ({}) NnonSat: {}".format(detNum, detName, len(nonSatArr)), fontsize=12)
    for ax in axes.flat:
        unitStr = "(scale factor)" if dataType == "flat" else "(counts?)"
        ax.set(xlabel=("{} " + unitStr).format(dataType), ylabel="Number")
        ax.label_outer()

    plt.suptitle("{} source: {} expId: {}".format(dataType, fileUsedDir, expId), fontsize=14, y=1.02)


def plotDiffHists(expId, butler, expId2, butler2, dataType, histColor=None, detectorOrder=None,
                  doPercentDiff=False, fileUsedDir=None, fileUsedDir2=None):
    fig, axes = plt.subplots(3,3, figsize=(13, 13), constrained_layout=True)
    nBins = 1000000 if (doPercentDiff and dataType == "bias") else 1000
    ymax = 5000000 if (doPercentDiff and dataType == "bias") else 200000
    if dataType == "flat" and doPercentDiff:
        nBins, ymax = 2000, 500000
    if dataType == "dark" and doPercentDiff:
        nBins, ymax = 500000, 20000000
    camera = butler.get("camera")
    histKwargs = dict(histtype='stepfilled', alpha=0.8, bins=nBins, fc=histColor, ec="k")
    print("det      min         max      median    mean    stddev      var")
    formatStr = "{:2d} {:11.2f} {:11.2f} {:8.2f} {:8.2f} {:8.2f} {:11.2f}"
    dataArr1DList = []
    xmin, xmax = 1e12, -1e12
    for iCcd, ax in enumerate(axes.flatten()):
        detNum = detectorOrder[iCcd]
        dataId = {'expId': expId, 'detector': detNum}
        dataImage = butler.get(dataType, dataId=dataId)
        dataArr = dataImage.maskedImage.image.array
        dataId2 = {'expId': expId2, 'detector': detNum}
        dataImage2 = butler2.get(dataType, dataId=dataId2)
        dataArr2 = dataImage2.maskedImage.image.array
        if doPercentDiff:
            diffArr = 200*((dataArr - dataArr2)/(dataArr + dataArr2))
        else:
            diffArr = dataArr - dataArr2
        diffArr = np.nan_to_num(diffArr, posinf=0.0, neginf=0.0)
        diffArr1D = diffArr.ravel()
        dataArr1DList.append(diffArr1D)
        percentile = 5 if doPercentDiff else 0.1
        xmin = min(xmin, np.percentile(diffArr1D, percentile))
        xmax = max(xmax, np.percentile(diffArr1D, 100 - percentile))
    for iCcd, ax in enumerate(axes.flatten()):
        detNum = detectorOrder[iCcd]
        detName = camera[detNum].getName()
        dataArr1D = dataArr1DList[iCcd]
        nanmin, nanmax = np.nanmin(dataArr1D), np.nanmax(dataArr1D)
        datamin, datamax = np.min(dataArr1D), np.max(dataArr1D)
        nanmedian, nanmean = np.nanmedian(dataArr1D), np.nanmean(dataArr1D)
        nanstd, nanvar = np.nanstd(dataArr1D), np.nanvar(dataArr1D)
        if iCcd == 0 and nanstd < 1.0:
            formatStr = "{:2d} {:9.3f} {:10.3f} {:8.4f} {:8.4f} {:8.4f} {:9.4f}"
        print(formatStr.format(iCcd, datamin, datamax, nanmedian, nanmean, nanstd, nanvar))
        ax.hist(dataArr1D, label="Det: {} ({})".format(detNum, detName), **histKwargs)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(0, ymax)
        ax.set_title("Det: {} ({})".format(detNum, detName, fontsize=12))
    percentStr = "Percent " if doPercentDiff else "Absolute "
    for ax in axes.flat:
        ax.set(xlabel=percentStr + "Difference: {}".format(dataType), ylabel="Number")
        ax.label_outer()
    plt.suptitle("{}: {}Difference Distributions ({} - {})".
                 format(dataType, percentStr, fileUsedDir, fileUsedDir2), fontsize=14, y=1.01)


def plotShowCamera(butler, expId, dataType, dataType2=None, butler2=None, expId2=None,
                   dataType3=None, butler3=None, expId3=None, binSize=10,
                   cmap="viridis", doVariance=False, doSignalToNoise=False, doDiffIm=False,
                   doPercentDiffIm=False, doAddIm=False, doAddThenSubtractIm=False,
                   fileUsedDir=None, fileUsedDir2=None, figureSideSize=12):
    subTypeStr = ""
    camera = butler.get("camera")
    plt.figure(figsize=(figureSideSize, figureSideSize))
    disp = afwDisplay.Display(1, "matplotlib")
    disp.scale("linear", "zscale")
    disp.setImageColormap(cmap)
    callback = None
    if doVariance:
        subTypeStr = " variance"
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True,
                                                  callback=lambda im, ccd, imageSource:
                                                  varianceOrSignalToNoiseCallback(im, ccd, imageSource,
                                                                                  dataType=dataType,
                                                                                  butler=butler, expId=expId,
                                                                                  doVariance=True))
    elif doSignalToNoise:
        subTypeStr = " S/N"
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True,
                                                  callback=lambda im, ccd, imageSource:
                                                  varianceOrSignalToNoiseCallback(im, ccd, imageSource,
                                                                                  dataType=dataType,
                                                                                  butler=butler, expId=expId,
                                                                                  doSignalToNoise=True))
    elif doDiffIm:
        subTypeStr = " Absolute Difference Image"
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True,
                                                  callback=lambda im, ccd, imageSource:
                                                  diffImCallback(im, ccd, imageSource, dataType=dataType,
                                                                 butler=butler, expId=expId,
                                                                 butler2=butler2, expId2=expId2,
                                                  doAbsoluteDiff=doDiffIm))

    elif doPercentDiffIm:
        subTypeStr = " Percent Difference Image"
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True,
                                                  callback=lambda im, ccd, imageSource:
                                                  diffImCallback(im, ccd, imageSource, dataType=dataType,
                                                                 butler=butler, expId=expId,
                                                                 butler2=butler2, expId2=expId2,
                                                                 doPercentDiff=doPercentDiffIm))
    elif doAddIm:
        subTypeStr = " Addition of two Images"
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True,
                                                  callback=lambda im, ccd, imageSource:
                                                  addImCallback(im, ccd, imageSource, dataType=dataType,
                                                                butler=butler, expId=expId,
                                                                dataType2=dataType2,
                                                                butler2=butler2, expId2=expId2))
    elif doAddThenSubtractIm:
        subTypeStr = " Addition of two Images and Subtraction of a Third"
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True,
                                                  callback=lambda im, ccd, imageSource:
                                                  addThenSubtractImCallback(im, ccd, imageSource,
                                                                            dataType=dataType,
                                                                            butler=butler, expId=expId,
                                                                            dataType2=dataType2,
                                                                            butler2=butler2, expId2=expId2,
                                                                            dataType3=dataType3,
                                                                            butler3=butler3, expId3=expId3))

    else:
        imageSource = cameraGeomUtils.ButlerImage(butler, dataType, expId=expId, verbose=True)
    if "Difference" in subTypeStr:
        titleStr = ("{}{} source1: {} source2: {}  BinSize={}".
                    format(dataType, subTypeStr, fileUsedDir, fileUsedDir2, binSize))
    else:
        titleStr = ("{}{} source: {}  expId: {}  BinSize={}".
                    format(dataType, subTypeStr, fileUsedDir, expId, binSize))
    mos = cameraGeomUtils.showCamera(camera, imageSource=imageSource, detectorNameList=None,
                                     binSize=binSize, display=disp, title=titleStr, ctype=afwDisplay.GREEN,
                                     textSize=3)


def varianceOrSignalToNoiseCallback(im, ccd, imageSource, dataType=None, butler=None, expId=None,
                                    doVariance=False, doSignalToNoise=False):
    dataId = {'expId': expId, 'detector': ccd.getId()}
    dataImage = butler.get(dataType, dataId=dataId)
    if doVariance:
         im[:] = dataImage.variance
    if doSignalToNoise:
         err = np.sqrt(dataImage.variance.array)
         signalToNoise = dataImage.image.array/err
         im.array = signalToNoise
    return im


def diffImCallback(im, ccd, imageSource, dataType=None, butler=None, expId=None, butler2=None, expId2=None,
                   doAbsoluteDiff=False, doPercentDiff=False):
    dataId = {'expId': expId, 'detector': ccd.getId()}
    dataImage = butler.get(dataType, dataId=dataId)
    dataId2 = {'expId': expId2, 'detector': ccd.getId()}
    dataImage2 = butler2.get(dataType, dataId=dataId2)
    if doAbsoluteDiff:
        diffArray = dataImage.image.array - dataImage2.image.array
        im.array = diffArray
    if doPercentDiff:
        diffArray = (200.0*(dataImage.image.array - dataImage2.image.array)
                     /(dataImage.image.array + dataImage2.image.array))
        diffArray = np.nan_to_num(diffArray, posinf=0.0, neginf=0.0)
        im.array = diffArray
    return im


def addImCallback(im, ccd, imageSource, dataType=None, butler=None, expId=None, dataType2=None,
                  butler2=None, expId2=None):
    dataId = {'expId': expId, 'detector': ccd.getId()}
    dataImage = butler.get(dataType, dataId=dataId)
    dataId2 = {'expId': expId2, 'detector': ccd.getId()}
    dataImage2 = butler2.get(dataType2, dataId=dataId2)
    diffArray = dataImage.image.array + dataImage2.image.array
    im.array = diffArray
    return im


def addThenSubtractImCallback(im, ccd, imageSource, dataType=None, butler=None, expId=None,
                              dataType2=None, butler2=None, expId2=None,
                              dataType3=None, butler3=None, expId3=None):
    dataId = {'expId': expId, 'detector': ccd.getId()}
    dataImage = butler.get(dataType, dataId=dataId)
    dataId2 = {'expId': expId2, 'detector': ccd.getId()}
    dataImage2 = butler2.get(dataType2, dataId=dataId2)
    dataId3 = {'expId': expId3, 'detector': ccd.getId()}
    dataImage3 = butler3.get(dataType3, dataId=dataId3)
    diffArray = dataImage.image.array + dataImage2.image.array - dataImage3.image.array
    im.array = diffArray
    return im
