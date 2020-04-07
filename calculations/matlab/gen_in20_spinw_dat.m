function hf = gen_in20_spinw_dat()

yn = 5.2; % y-min
ym = 30;  % y-max
load 'in20_fitted_peaks.mat';
lns = [];
for ist=1:length(setindex{3})
    hkl=sscanf(setindex{3}{ist},'(%f %f %f)');
    ihkl=sum(abs(cumhkl(:,1:3)-(repmat(hkl',length(cumhkl),1)))'); 
    ihkl=find(ihkl==min(ihkl));
    if(sum(abs(cumhkl(ihkl(1),1:3)-hkl'))<1e-1)
        for iu=unique(cumhkl(ihkl,4))'
            ln = [iu hkl' yn ym];
            for ie=1:(length(fits{ist}{1})-1)/3
                if(abs(fits{ist}{1}((ie-1)*3+1))<1e7 && fits{ist}{1}((ie-1)*3+2)>yn && fits{ist}{1}((ie-1)*3+2)<ym)
                    if (iu == 2.991 && fits{ist}{1}((ie-1)*3+2) > 20)
                        continue; % Remove a spurious peak
                    end
                    area = abs(fits{ist}{1}((ie-1)*3+1)) / 2e5;
                    centre = fits{ist}{1}((ie-1)*3+2);
                    sigma = fits{ist}{1}((ie-1)*3+3) / 2.35;
                    ln = [ln area centre sigma]; 
                end
            end
            %if iu == 0.4185  % Has 3 modes - SpinW fitspec fits only 2 modes
            %    ln = [ln(1:6) mean(ln([7 10])) mean(ln([8 11])) sum(ln([9 12])) ln(13:15)];
            %end
            if numel(ln) > 6
                if numel(ln) < 15
                    ln((end+1):15) = 0;
                end
                %disp(ln);
                lns = [lns; ln];
            end
        end
    end
end

[~, ids] = sort(lns(:,1));
lns(ids, 1:end)

fid = fopen('in20_fitted_modes.txt', 'w')
fprintf(fid, '   QH        QK        QL        EMIN     EMAX       I1       EN1        sigma1    I2       EN2        sigma2    I3        EN3       sigma3\n');
%fprintf(fid, '   QH        QK        QL        EMIN     EMAX       I1       EN1        sigma1    I2       EN2        sigma2\n');
fprintf(fid, [repmat('%9.4f ', 1, 14) '\n'], lns(ids, 2:end)');
fclose(fid);
