function solCuSeq( matLstCuCx )  
    for it = 0 : div(20000, itLm)
        matCu = cu(matLstCatCx[:,:,it*itLm+1:min((it+1)*itLm,20000)]);
        solCuCx = CUDA.CUSOLVER.heevjBatched!( 'V','U', matCu );
    end                              
end


function solCuSeqBlock( matLstCuCx )  
    for it = 0 : div(20000, itLm)
        solCuCx = CUDA.CUSOLVER.heevjBatched!( 'V','U', matLstCuCx[:,:,it*itLm+1:min((it+1)*itLm,20000)] );
    end                              
end


function solCuSeqView( matLstCuCx )  
    for it = 0 : div(20000, itLm)
        solCuCx = CUDA.CUSOLVER.heevjBatched!( 'V','U', view( matLstCuCx, :,:,it*itLm+1:min((it+1)*itLm,20000) ) );
    end                              
end


@time divB_profile( 3, [20],1000, [80,16,16], 1000, fileNameMod = "threadedTest" );