      real function round ( x )
!  from  * a practical guide to splines *  by c. de boor
!alled in example 1  of chapter xiii
      real x,   flip,size
      common /rount/ size
      data flip /-1./
      flip = -flip
      round = x + flip*size
                                        return
      end
