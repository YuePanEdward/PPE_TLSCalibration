function plot_scanner(scanner_ep,scanner_id,fig_index)
%plot_scanner: draw a scanner's pose
%inputs:
% - scanner_ep: 6dof ep of the scanner
% - scanner_id: id of the scanner
% - fig_index is the figure number

axis_length = 1.0; 
offset = 2.5;
linewidth = 3;

po=scanner_ep(4:6);
R_so=eul2rotm(scanner_ep(1:3)','XYZ');

po_axis = R_so' * (axis_length.* eye(3)) + repmat(po,[1,3]);
px = po_axis(:,1);
py = po_axis(:,2);
pz = po_axis(:,3);

figure(fig_index);
plot3(po(1),po(2),po(3),'k.','markersize',30);

hold on, text(po(1)+offset*(0.5-rand(1)), po(2)+offset*(0.5-rand(1)), po(3)+offset*(0.5-rand(1)), ['S',num2str(scanner_id)], 'Color', 'k','FontSize', 14);
hold on, line([po(1); px(1)], [po(2); px(2)], [po(3); px(3)], 'Color', 'r', 'LineWidth', linewidth, 'DisplayName','Scanner_axis_x'); 
hold on, line([po(1); py(1)], [po(2); py(2)], [po(3); py(3)], 'Color', 'g', 'LineWidth', linewidth, 'DisplayName','Scanner_axis_y'); 
hold on, line([po(1); pz(1)], [po(2); pz(2)], [po(3); pz(3)], 'Color', 'b', 'LineWidth', linewidth, 'DisplayName','Scanner_axis_z'); 
axis equal;
grid on;

end