clearvars;

%length of input
k = input('Length of input to encoder');

%number of times a bit should be repeated
n=3:2:35;

%rate
rate=1./n;

%generate input
in = gen_input(k);

%sigma
sigma=1;

%number of times to decode
Nsim=10000;


for j=1:17

    %set initial error of each channel to 0
    error_gaussian=0;

    for i=1:Nsim

        %generate generator matrix
        generator_matrix=gen_generator_matrix(k,n(j));

        %encode input
        encoded_input=encode_input(generator_matrix,in);

        %pass encoded bits through channels
        channel_pass_gaussian = gaussian_out(encoded_input,sigma);

        %decode
        output_gaussian=decode(channel_pass_gaussian,k,n(j),"Gaussian");

        %calculate error
        error_gaussian=error_gaussian+compare(in,output_gaussian,k);

    end

    %calculate error rate
    error_rate_gaussian(j)=error_gaussian/(k*Nsim);

end

figure();
semilogy(n,error_rate_gaussian,'^-','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]);
xlabel('Number of times a bit is being repeated'); 
ylabel('Probability of Bit Error');  
grid on;

legend('Gaussian'); 
axis([0 35 1e-6 1]); 
set(gca,'xtick',3:2:35);
    





function [in] = gen_input(k)
    p = 0.5;
    a = rand(1,k);
    for i = 1:k
        a(i) = ceil(a(i)-p);
    end
    in = a;
end

function [generator_matrix] = gen_generator_matrix (k,n)
    gen_matrix=zeros(n*k,k);
    t=1;
    for j=1:k
        for i=t:t+n-1
            gen_matrix(i,j)=1;
        end
        t=t+n;
    end
    generator_matrix=gen_matrix;
end

function [encoded_input] = encode_input(generator_matrix,in)
    out=generator_matrix*in';
    encoded_input=out';
end

function [gaussian] = gaussian_out(encoded_input,sigma)
    for i = 1:size(encoded_input,2)
        n = randn*sigma;
        encoded_input(i) = (2*encoded_input(i)-1) + n;
    end
    gaussian = encoded_input;
end

function [output] = decode(channel_out,k,n,channel)
    out=zeros(1,k);
    if(channel=="Gaussian")
        for i=1:k
            eucledian_dist_1=0;
            eucledian_dist_neg_1=0;
            for j=(i-1)*n+1:i*n
                eucledian_dist_1=eucledian_dist_1+(channel_out(j)-1)*(channel_out(j)-1);
                eucledian_dist_neg_1=eucledian_dist_neg_1+(channel_out(j)-(-1))*(channel_out(j)-(-1));
            end
            if(eucledian_dist_1>eucledian_dist_neg_1)
                out(i)=0;
            else
                out(i)=1;
            end
        end
    end
    output=out;
end

function [error] = compare(in,output,k)
    err=0;
    for i=1:k
        if(in(i)~=output(i))
            err=err+1;
        end
    end
    error=err;
end