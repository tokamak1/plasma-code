WriterObj=VideoWriter(['Lyapunov exponent t = 0 ~ ',num2str(timep(end)),...
                       'dt = ',num2str(timep(1)),'.avi']);
WriterObj.FrameRate = 5;
open(WriterObj);
writeVideo(WriterObj,F3);
close(WriterObj);