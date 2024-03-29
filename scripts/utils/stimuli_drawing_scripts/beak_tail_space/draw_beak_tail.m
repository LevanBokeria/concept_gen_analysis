function [rBeak, rBody, rTail, rBeakTip, rTailTip] ...
    = draw_beak_tail(SF_all,SF_body,...
               xAdj,yAdj,x,y)

% Re-generate the parameters
[beak_height, beakTip_width, beakTip_height, ...
   tail_height, tailTip_width, ...
   bird_height, bird_width, ...
   beak_body_offset_y, beak_body_offset_x, ...
   tail_body_offset, beakTip_offset, tailTip_offset, ...
   overlap] ...
   = parameters_beak_tail(SF_all,SF_body);

% Define the boxes
rBody = [0, ...
   0, ...
   bird_width, ...
   bird_height];

rBody = CenterRectOnPoint(rBody, x,y);

rBeak = [rBody(1) + overlap - xAdj - beak_body_offset_x, ...
   rBody(2) + beak_body_offset_y, ...
   rBody(1) + overlap - beak_body_offset_x, ...
   rBody(2) + beak_body_offset_y + beak_height];

rBeakTip = [rBeak(1) - beakTip_width + overlap, ...
   rBeak(2) + beakTip_offset, ...
   rBeak(1) + overlap, ...
   rBeak(2) + beakTip_offset + beakTip_height];

rTail = [rBody(3), ...
   rBody(2) + tail_body_offset, ...
   rBody(3) + yAdj, ...
   rBody(2) + tail_body_offset + tail_height];

rTailTip = [rTail(3) - overlap, ...
   rTail(4) - beakTip_height + tailTip_offset, ...
   rTail(3) - overlap + tailTip_width, ...
   rTail(4) + tailTip_offset];

end