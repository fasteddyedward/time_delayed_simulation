function save_movie(MovieVector,movie_name,frame_rate)
myWriter=VideoWriter(movie_name);
myWriter.FrameRate=frame_rate;
open(myWriter);
writeVideo(myWriter,MovieVector);
close(myWriter);
end