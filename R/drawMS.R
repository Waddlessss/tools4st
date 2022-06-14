
#' drawMS
#' @description Draw an averaged MS scan during a certain retention time.
#'
#' @param rawDataName character. File name of the raw data (mzML or mzXML).
#' @param mzRange vector. m/z range to plot.
#' @param rtRange vector. Retention time to average.
#' @param mzTol float. m/z tolerance to average spectra.
#' @param msLevel int. MS level (1 or 2).
#' @param mzLabelDis float. A value to specify the range that only one m/z will be labeled, in Da.
#' @param intLabelTol float. Relative intensity larger than this value can be labeled.
#' @param dcm int. Number of decimals for m/z labels.
#' @param res int. Resolution of the output figure.
#'
#' @return
#' This function outputs a figure in the current working directory.
#'
#' @import MSnbase
#'
#' @export
#'
#' @examples
#' # drawMS(rawDataName)


drawMS = function(rawDataName, mzRange=c(700,850), rtRange=c(0,Inf), mzTol=0.005, msLevel=1,
                  mzLabelDis=20, intLabelTol=3, res=600, dcm=2){

  # Read raw MS data
  rawData = readMSData(files=rawDataName, msLevel.=msLevel)

  # Filter MS scans by defined retention time window
  rtFilter = rtime(rawData) <= rtRange[2]*60 & rtime(rawData) >= rtRange[1]*60
  mzList = mz(rawData)[rtFilter]
  intList = intensity(rawData)[rtFilter]

  # Generate vectors to store the averaged scans
  mzPlot = c(-Inf)
  intPlot = c(-Inf)
  counter = c(1)

  # Check if mzList is empty
  if (length(mzList) == 0) {
    message("No MS scan detected in defined retention time window.")
    return(NA)
  }

  # Average the selected scans
  for (i in 1:length(mzList)) {
    mzFilter = mzList[[i]] <= mzRange[2] & mzList[[i]] >= mzRange[1]
    mz = mzList[[i]][mzFilter]
    intensity = intList[[i]][mzFilter]
    # intensity = intensity / max(intensity) * 100

    for (j in 1:length(mz)) {
      temp = abs(mz[j]-mzPlot)
      if (any(temp < mzTol)) {
        idx = which.min(temp)
        mzPlot[idx] = (mzPlot[idx]*counter[idx] + mz[j])/(counter[idx]+1)
        intPlot[idx] = (intPlot[idx]*counter[idx] + intensity[j])/(counter[idx]+1)
        counter = counter + 1
      } else {
        mzPlot = c(mzPlot, mz[j])
        intPlot = c(intPlot, intensity[j])
        counter = c(counter, 1)
      }
    }
  }
  mzPlot = mzPlot[-1]
  intPlot = intPlot[-1]/max(intPlot[-1])*100
  # Plot
  windowsFonts(
    A=windowsFont("Arial")
  )

  png("image.png", res=600, width = 4.5, height = 4.5, units = 'in')
  par(family="A", bg="transparent")
  plot(mzPlot, intPlot, type="h", xlab = expression(italic("m/z")), ylab="R.A.",
       cex.lab = 1.5, cex.axis = 1.2, las=1, ylim=c(0,110))
  mzToLabel = findLabelText(mzPlot, intPlot, mzLabelDis = mzLabelDis, intTol = intLabelTol)
  text(mzToLabel[[1]], mzToLabel[[2]], round(mzToLabel[[1]],dcm), pos = 3)
  dev.off()
}


findLabelText = function(mzPlot, intPlot, labNumber=4, intTol=3, mzLabelDis=20){
  mzLabel = c()
  intLabel = c()
  while (max(intPlot) > intTol & length(mzLabel)<labNumber) {
    mzLabel = c(mzLabel, mzPlot[which.max(intPlot)])
    intLabel = c(intLabel, intPlot[which.max(intPlot)])
    intPlot[abs(mzPlot-mzPlot[which.max(intPlot)])<mzLabelDis] = 0
  }
  return(list(mzLabel, intLabel))
}





