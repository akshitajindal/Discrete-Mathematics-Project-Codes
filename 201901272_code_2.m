clearvars;

%length of input
k = input('Length of input to encoder');

%number of times a bit should be repeated
n=3:2:35;

%rate
rate=1./n;

%generate input
in = gen_input(k);

%probability
p=0.1;

%number of times to decode
Nsim=10000;

for j=1:17
    p_error_bsc(j)=0;
    for t=(n(j)+1)/2:n(j)
        error=nchoosek(n(j),t)*power(p,t)*power((1-p),(n(j)-t));
        p_error_bsc(j) = p_error_bsc(j) + error;
    end
    p_error_bec(j)=p^n(j);
end


for j=1:17

    %set initial error of each channel to 0
    error_bsc=0;
    error_bec=0;
    error_gaussian=0;

    for i=1:Nsim

        %generate generator matrix
        generator_matrix=gen_generator_matrix(k,n(j));

        %encode input
        encoded_input=encode_input(generator_matrix,in);

        %pass encoded bits through channels
        channel_pass_bsc = bsc_out(encoded_input,p);
        channel_pass_bec = bec_out(encoded_input,p);

        %decode
        output_bsc=decode(channel_pass_bsc,k,n(j),"BSC");
        output_bec=decode(channel_pass_bec,k,n(j),"BEC");

        %calculate error
        error_bsc=error_bsc+compare(in,output_bsc,k);
        error_bec=error_bec+compare(in,output_bec,k);

    end

    %calculate error rate
    error_rate_bsc(j)=error_bsc/(k*Nsim);
    error_rate_bec(j)=error_bec/(k*Nsim);

end

figure();
semilogy(n,error_rate_bsc,'o-','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
semilogy(n,p_error_bsc,'o:','linewidth',2,'markerfacecolor','b','markeredgecolor','b');
hold on;
xlabel('Number of times a bit is being repeated'); 
ylabel('Probability of Bit Error'); 

semilogy(n,error_rate_bec,'d-','linewidth',2,'color',[0 0.4 0.9],'markerfacecolor',[0 0.4 0.9],'markeredgecolor',[0 0.4 0.9]);
hold on;
semilogy(n,p_error_bec,'^:','linewidth',2,'color',[0 0.5 0],'markerfacecolor',[0 0.5 0],'markeredgecolor',[0 0.5 0]); 
hold on;
xlabel('Number of times a bit is being repeated'); 
ylabel('Probability of Bit Error'); 
grid on;

legend('BSC Experimental','BSC Theoritical','BEC Experimental','BEC Theoritical'); 
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

function [bsc] = bsc_out(encoded_input,p)
    for i = 1:size(encoded_input,2)
        x = rand;
        if(x<p)
            if(encoded_input(i)==0)
                encoded_input(i)=1;
            elseif(encoded_input(i)==1)
                encoded_input(i)=0;
            end
        end
    end
    bsc = encoded_input;
end

function [bec] = bec_out(encoded_input,p)
    for i = 1:size(encoded_input,2)
        x = rand;
        if(x<p)
            %let the erased bit be equal to 5
            encoded_input(i) = 5;
        end
    end
    bec = encoded_input;
end

function [output] = decode(channel_out,k,n,channel)
    out=zeros(1,k);
    if(channel=="BSC")
        for i=1:k
            no_zeros=0;
            no_ones=0;
            for j=(i-1)*n+1:i*n
                if(channel_out(j)==0)
                    no_zeros=no_zeros+1;
                elseif(channel_out(j)==1)
                    no_ones=no_ones+1;
                end
            end
            if(no_zeros>no_ones)
                out(i)=0;
            else
                out(i)=1;
            end
        end
    end
    if(channel=="BEC")
        for i=1:k
            no_zeros=0;
            no_ones=0;
            erased_bits=0;
            for j=(i-1)*n+1:i*n
                if(channel_out(j)==0)
                    no_zeros=no_zeros+1;
                elseif(channel_out(j)==1)
                    no_ones=no_ones+1;
                elseif(channel_out(j)==5)
                    erased_bits=erased_bits+1;
                end
            end
            if(no_zeros==0 && no_ones~=0)
                out(i)=1;
            elseif(no_ones==0 && no_zeros~=0)
                out(i)=0;
            elseif(no_zeros==0 && no_ones==0)
                out(i)=5; %all bits got erased therefore output bit will be erased
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