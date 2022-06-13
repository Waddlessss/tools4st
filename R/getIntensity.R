
#' getIntensity
#'
#' @description Extract intensities by m/z value and retention time
#'
#' @param rawDataName character. File name of raw MS data (mzML or mzXML).
#' @param ionTableName character. A csv file specifying the m/z values and retention time windows.
#' @param mzTol float. m/z tolerance to average spectra.
#' @param msLevel int. MS level (1 or 2).
#' @param output Boolean. Set to \code{TRUE} then output the result table.
#'
#' @return
#' This function either returns the result table or output it to csv file.
#' @export
#'
#' @examples
#' # getInensity(rawDataName, ionTableName)

getIntensity = function(rawDataName, ionTableName, mzTol=0.005, msLevel=1, output=TRUE){
  rawData = MSnbase::readMSData(files=rawDataName, msLevel.=msLevel)
  ionTable = read.csv(ionTableName)

  mz = as.numeric(ionTable[,1])
  rt = ionTable[,2]
  mz = mz[mz!=""]
  rt = rt[rt!=""]

  resultTable = data.frame(matrix(nrow = length(mz), ncol = length(rt)))
  rownames(resultTable) = round(mz,4)
  colnames(resultTable) = rt

  for (i in 1:length(rt)) {
    rtRange = as.numeric(unlist(strsplit(rt[i], split = ";")))
    rtFilter = rtime(rawData) <= rtRange[2]*60 & rtime(rawData) >= rtRange[1]*60
    mzList = mz(rawData)[rtFilter]
    intList = intensity(rawData)[rtFilter]

    for (j in 1:length(mz)) {
      intSeq = c()
      for (k in 1:length(mzList)) {
        mzDiff = abs(mz[j]-mzList[[k]])
        idx = which.min(mzDiff)
        if (mzDiff[idx] < mzTol) {
          intSeq = c(intSeq, intList[[k]][idx])
        }
      }
      resultTable[j,i] = mean(intSeq)
    }
  }

  if (output) {
    outputFileName = paste0(tools::file_path_sans_ext(rawDataName), "_extracted.csv")
    write.csv(resultTable, outputFileName)
  } else {
    return(resultTable)
  }
}
