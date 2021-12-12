tmap = imread('IMG_20201007_125027-WangCenterRight.jpg');
R=tmap(:,:,1);
G=tmap(:,:,2);
B=tmap(:,:,3);

%imshow(R);
%figure,imshow(G);
%figure,imshow(B);

[r,c]=size(R);
I=uint8(zeros(r,c));
for i=1:r
    for j=1:c
        gv=0.299 * R(i,j) + 0.587 * G(i,j) + 0.114 * B(i,j);
        %fprintf("%d\n",gv);
        I(i,j)=gv;
    end
end

%figure
%imshow(I);  %gray level image I
%title('gray level');

P=zeros(255,1);  %h(k)
for i=1:r
    for j=1:c
        
        ind=I(i,j);
        %fprintf("index is %d\n",ind);
        temp=P(ind+1,1);
        P(ind+1,1)=temp+1;
    end
end

imgsize=double(r*c);

probdev=P./imgsize;
H=uint8(zeros(r,c));  %equalized matrix

for i=1:r
    for j=1:c
        prob=I(i,j);
        result=0;
        for k=1:prob
            result=result+probdev(k);
        end
        result=uint8(result*255);
        H(i,j)=result;
    end
    
end
figure, imshow(H),title('equalized'); 

formatSpec = '%d %d';
sizeA = [5 5];
fileID1 = fopen('filter1.txt','r');
filter1 = fscanf(fileID1,formatSpec,sizeA);
%disp(filter1);

formatSpec = '%f %f';
fileID2 = fopen('filter2.txt','r');
filter2 = fscanf(fileID2,formatSpec,sizeA);
%disp(filter2);

%convolution
C_1=uint8(zeros(r,c));

for i=3:r-2
    for j=3:c-2
        temp=double(0);
        for k=1:5
           for l=1:5
               v=filter1(k,l)*double(I(i-(k-3),j-(l-3)));
               temp=temp+v;
                
           end
        end
        %fprintf('temp is %d\n',temp);
        C_1(i,j)=uint8(temp);
    end
    
end

min1=0;
max1=0;

for i=1:r   %find min and max
    for j=1:c
        if C_1(i,j)>max1
            max1=C_1(i,j);
        elseif C_1(i,j)<min1
            min1=C_1(i,j);
        end
        
    end
end
dif=max1-min1;
for i=1:r  %reset image intensity
    for j=1:c
        int_min=double(C_1(i,j)-min1);
        inttemp=double(int_min*255);
        %fprintf('intensity-min is %d\n',int_min);
        %fprintf('intensity is %d\n',inttemp);
        tempvalue=uint8(inttemp/dif);
        %fprintf('%d\n',tempvalue);
        %C_1(i,j)=tempvalue;
        
    end
end

figure, imshow(C_1),title('filter1');

C2=uint8(zeros(r,c));

for i=3:r-2
    for j=3:c-2
        temp=0;
        for k=1:5
           for l=1:5
               temp=temp+filter2(k,l)*I(i-(k-3),j-(l-3));
           end
        end
        C2(i,j)=temp;
    end
    
end

min2=0;
max2=0;

for i=1:r   %find min and max
    for j=1:c
        if C2(i,j)>max2
            max2=C2(i,j);
        elseif C2(i,j)<min2
            min2=C2(i,j);
        end
        
    end
end
dif2=max2-min2;
for i=1:r  %reset image intensity
    for j=1:c
        int_min=uint32(C2(i,j)-min2);
        inttemp=uint32(uint32(int_min)*255);
        %fprintf('intensity-min is %d\n',int_min);
        %fprintf('intensity is %d\n',inttemp);
        tempvalue=uint8(double(inttemp)/double(dif2));
        %fprintf('%d\n',tempvalue);
        C2(i,j)=tempvalue;
        
    end
end
figure, imshow(C2),title('filter2');



%gaussian filter

sigma=2.0;
summ=0; 
g=zeros(9,1);
for k = 1 : 9
    g(k,1) =  exp( -((k-5)*(k-5)) / (2*sigma*sigma) );
    summ = summ+g(k,1);
 end
for k = 1 : 9
   g(k,1) =  g(k,1)/summ;
end



h1=uint8(zeros(r,c));
h2=uint8(zeros(r,c));

for i=1:r   %find min and max
    for j=1:c
        h1(i,j)=I(i,j);
        h2(i,j)=I(i,j);
    end
end

for i = 1: r 
    for j = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 9 
           summ = summ + g(k,1)* I(i,j-(k-5)); 
        end 
        h1(i,j) = summ; 
     end 
end 

%filter each column j
for j = 1: r 
    for i = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 9 
           summ = summ + g(k,1)* h1(i-(k-5),j); 
        end 
        h2(i,j) = summ; 
    end 
end

figure, imshow(h2),title('gaussian');  %matrix after gaussian filtering

ts=5;  %threshold
edge=uint8(zeros(r,c));
for i = 2: (r-1) 
    for j = 2 : (c-1) 
        dx=double(h2(i+1,j))-double(h2(i,j));
        dy=double(h2(i,j+1))-double(h2(i,j));
        dir=double(dx^2+dy^2);
        mag=sqrt(dir);
        %fprintf('%d\n',mag);
        if mag>ts
            edge(i,j)=255;
        end
    end 
end

figure, imshow(edge) ,title('edge');


%corner detection

A=double(zeros(r,c));
B=double(zeros(r,c));
C=double(zeros(r,c));

for i=1:(r-1)   %calulate A B C
    for j=1:(c-1)
        dx=double(h2(i+1,j))-double(h2(i,j));
        dy=double(h2(i,j+1))-double(h2(i,j));
        ctemp=dx*dy;
        A(i,j)=dx;
        B(i,j)=dy;
        C(i,j)=ctemp;
        
    end
end
figure, imshow(A),title('A');
figure, imshow(B),title('B');
figure, imshow(C),title('C');
%////////////////////////////////////////////////////////////////////////////
sigma=5.5;
summ=0; 
g=zeros(11,1);
for k = 1 : 11
    g(k,1) =  exp( -((k-6)*(k-6)) / (2*sigma*sigma) );
    summ = summ+g(k,1);
 end
for k = 1 : 11
   g(k,1) =  g(k,1)/summ;
end

A1=double(zeros(r,c));
A2=double(zeros(r,c));

for i=1:r   
    for j=1:c
        A1(i,j)=A(i,j);
        A2(i,j)=A(i,j);
    end
end

for i = 1: r 
    for j = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 11 
           summ = summ + g(k,1)* A(i,j-(k-6)); 
        end 
        A1(i,j) = summ; 
     end 
end 

%filter each column j
for j = 1: r 
    for i = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 11 
           summ = summ + g(k,1)* A1(i-(k-6),j); 
        end 
        A2(i,j) = summ; 
    end 
end

figure, imshow(A2), title('A/');
%/////////////////////////////////////////////////////////////////
%////////////////////////////////////////////////////////////////////////////
sigma=5.5;
summ=0; 
g=zeros(11,1);
for k = 1 : 11
    g(k,1) =  exp( -((k-6)*(k-6)) / (2*sigma*sigma) );
    summ = summ+g(k,1);
 end
for k = 1 : 11
   g(k,1) =  g(k,1)/summ;
end

B1=double(zeros(r,c));
B2=double(zeros(r,c));

for i=1:r   %find min and max
    for j=1:c
       B1(i,j)=B(i,j);
        B2(i,j)=B(i,j);
    end
end

for i = 1: r 
    for j = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 11 
           summ = summ + g(k,1)* B(i,j-(k-6)); 
        end 
        B1(i,j) = summ; 
     end 
end 

%filter each column j
for j = 1: r 
    for i = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 11 
           summ = summ + g(k,1)* B1(i-(k-6),j); 
        end 
        B2(i,j) = summ; 
    end 
end

figure, imshow(B2), title('B/');
%/////////////////////////////////////////////////////////////////
%////////////////////////////////////////////////////////////////////////////
sigma=5.5;
summ=0; 
g=zeros(11,1);
for k = 1 : 11
    g(k,1) =  exp( -((k-6)*(k-6)) / (2*sigma*sigma) );
    summ = summ+g(k,1);
 end
for k = 1 : 11
   g(k,1) =  g(k,1)/summ;
end
disp(g);
C1=double(zeros(r,c));
C2=double(zeros(r,c));

for i=1:r   %find min and max
    for j=1:c
        C1(i,j)=C(i,j);
        C2(i,j)=C(i,j);
    end
end

for i = 1: r 
    for j = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 11 
           summ = summ + g(k,1)* C(i,j-(k-6)); 
        end 
        C1(i,j) = summ; 
     end 
end 

%filter each column j
for j = 1: r 
    for i = 6 : (c-5) 
        summ = 0; 
        for k = 1 : 11 
           summ = summ + g(k,1)* C1(i-(k-6),j); 
        end 
        C2(i,j) = summ; 
    end 
end

figure, imshow(C2), title('C/');
%/////////////////////////////////////////////////////////////////
RR=double(zeros(r,c));
corner=uint8(zeros(r,c));
for i=1:r 
    for j=1:c
        M=[A2(i,j) C2(i,j);
            C2(i,j) B2(i,j)];
        d=det(M);
        t=trace(M);
        RR(i,j)=d-0.04*(t^2);
        if RR(i,j)>0.05
            corner(i,j)=255;
        end
        %fprintf('%d\n',double(RR(i,j)));
    end
end

figure, imshow(corner), title('corner');

NM=uint8(zeros(r,c));
for i=2:(r-1)
    for j=2:(c-1)
        o=RR(i,j);
        up=RR(i,j-1);
        upleft=RR(i-1,j-1);
        upright=RR(i+1,j-1);
        left=RR(i-1,j);
        right=RR(i+1,j);
        down=RR(i,j+1);
        downleft=RR(i-1,j+1);
        downright=RR(i+1,j+1);
        %fprintf("%d %d %d %d %d %d %d %d\n",up,upleft,upright,left,right,downleft,down,downright)
        if (o>up)&&(o>down)&&(o>left)&&(o>right)&&(o>upright)&&(o>upleft)&&(o>downright)&&(o>downleft)
            NM(i,j)=255;
            %fprintf('success\n');
        end
    end
end

figure, imshow(NM), title('Non Max');

hn=zeros(1,8);
h=zeros(1,8);

for i=1:r
    for j=1:c
        if NM(i,j)==255
            fprintf('corner at %d, %d',i,j);
            dx=double(H(i+1,j))-double(H(i,j));
            dy=double(H(i,j+1))-double(H(i,j));
            %fprintf('dx dy are %d, %d',dx,dy);
            dir=atan2d(double(dy),double(dx));
            if dir<0
                dir=dir+360;
            end
            if (dir>=0)&&(dir<22.5)
                h(1)=h(1)+1;
            elseif(dir>=22.5)&&(dir<67.5)
                h(2)=h(2)+1;
            elseif(dir>=67.5)&&(dir<112.5)
                h(3)=h(3)+1;
            elseif(dir>=112.5)&&(dir<157.5)
                h(4)=h(4)+1;
            elseif(dir>=157.5)&&(dir<202.5)
                h(5)=h(5)+1;
            elseif(dir>=202.5)&&(dir<247.5)
                h(6)=h(6)+1;
            elseif(dir>=247.5)&&(dir<292.5)
                h(7)=h(7)+1;
            elseif(dir>=292.5)&&(dir<337.5)
                h(8)=h(8)+1;
            elseif(dir>=337.5)&&(dir<360)
                h(1)=h(1)+1;
            end
            disp(h);
            %fprintf('\ndir is %d\n',dir);
            
            
        end
    end
end

max=1;
for i=2:8
    if h(i)>h(max)
        max=i;
    end
end
%{
 h=h';
if max>5
    position=8-(max-5);
    h=circshift(h,position);
elseif max<5
    position=max-5;
    h=circshift(h,position);
end
h=h';
fprintf('Final histogram is: ');
disp(h);
%}
kk=5;
for i= 1 : 8
  hn(mod(kk-2+i,8)+1) = h(mod(max-2+i,8)+1);
end


fprintf('Final normalized histogram is: ');
disp(hn);
    














