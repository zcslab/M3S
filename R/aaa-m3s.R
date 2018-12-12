setClass('M3SMethod', representation = representation('VIRTUAL', modelFunction = 'function'))

setGeneric('m3s', function(x, method, ...) {standardGeneric('m3s')})

setMethod('m3s', c('matrix','M3SMethod'),
          function(x, method, ...) {
            MYCALL<-match.call()
            ret<-method@modelFunction(x,...)
            ret@Parameters<-c(list(Call=MYCALL, Method=method))
            return(ret)
          })

setMethod('m3s', c('matrix','function'),
          function(x,method, ...) {
            method <- method()
            biclust(x, method, ...)
          })

setMethod('m3s', c('matrix','character'),
          function(x,method, ...) {
            method <- get(method[1], mode="function")
            biclust(x,method, ...)
          })


setClass('M3S', representation = representation(Parameters = 'list', KS = 'list', info = 'list')
)

M3SResult <- function(mypara, ks, info) {
  return(new('M3S', Parameters = mypara, KS = ks, info = info))
}

setClass('P', contains = 'M3SMethod', prototype = prototype(modelFunction = function(x,minr=2,minc=2,number=100){bimaxbiclust(x,minr,minc,number)}))

P <- function() {
  return(new('P'))
}


