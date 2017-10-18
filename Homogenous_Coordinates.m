clear all
figure(1)
plot([1,1],[1,600])
hold on;
plot([1,800],[600,600])
hold on;
plot([800,800],[600,1])
hold on;
plot([800,800],[1,1])
e2p([1,2])
hold on;
[u,v] = ginput(4)
plot(u,v,'o')
hold on;
 
hom_c = e2p([u,v]')

line1 = cross(hom_c(:,1),hom_c(:,2))
line2 = cross(hom_c(:,3),hom_c(:,4))

line_l = cross([1,1,1]',[1,600,1]')
line_b = cross([1,600,1]',[800,600,1]')
line_r = cross([800,600,1]',[800,1,1]')
line_u = cross([1,1,1]',[800,1,1]')

l1_c(:,1) = p2e(cross(line1,line_l))
l1_c(:,2) = p2e(cross(line1,line_b))
l1_c(:,3) = p2e(cross(line1,line_r))
l1_c(:,4) = p2e(cross(line1,line_u))



pt = [[0;0],[0;0]]
w = 1
for i = 1:4
    disp(l1_c(:,i))
    a = check_point_inside(l1_c(:,i))
    if(a == 1)
        pt(:,w) = l1_c(:,i)
        w = w + 1
        if(w == 3)
            break
        end
    end
end    
plot(pt(1,:),pt(2,:))
hold on;


l2_c(:,1) = p2e(cross(line2,line_l))
l2_c(:,2) = p2e(cross(line2,line_b))
l2_c(:,3) = p2e(cross(line2,line_r))
l2_c(:,4) = p2e(cross(line2,line_u))



pt1 = [[0;0],[0;0]]
w = 1
for i = 1:4
    disp(l2_c(:,i))
    a = check_point_inside(l2_c(:,i))
    if(a == 1)
        pt1(:,w) = l2_c(:,i)
        w = w + 1
        if(w == 3)
            break
        end
    end
end    
plot(pt1(1,:),pt1(2,:))
hold on;
%Intersection of the two lines

int_pt = p2e(cross(line1,line2))
if(check_point_inside(int_pt) == 1)
    plot(int_pt(1),int_pt(2),'ro')
end

%Applying Homographies

figure(2)
h_lu = p2e(compute_homography(e2p([1,1]')));
h_lb = p2e(compute_homography(e2p([1,600]')));
h_rb = p2e(compute_homography(e2p([800,600]')));
h_ru = p2e(compute_homography(e2p([800,1]')));

plot([h_lu(1) h_lb(1)],[h_lu(2) h_lb(2)]);
hold on;
plot([h_lb(1) h_rb(1)],[h_lb(2) h_rb(2)]);
%plot([1,800],[600,600])
hold on;
plot([h_rb(1) h_ru(1)],[h_rb(2) h_ru(2)]);
%plot([800,800],[600,1])
hold on;
plot([h_ru(1) h_lu(1)],[h_ru(2) h_lu(2)]);
%plot([800,800],[1,1])
%e2p([1,2])
hold on;

a = p2e(compute_homography(e2p([u,v]')))

plot(a(1,:),a(2,:),'ro')
hold on;
pt = p2e(compute_homography(e2p(pt)))
pt1 = p2e(compute_homography(e2p(pt1)))

plot(pt(1,:),pt(2,:))
hold on;
plot(pt1(1,:),pt1(2,:))
hold on;

int_pt = p2e(compute_homography(e2p(int_pt)))
if(check_point_inside(int_pt) == 1)
    plot(int_pt(1),int_pt(2),'ro')
end