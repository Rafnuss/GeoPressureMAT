function [outputArg1,outputArg2] = classifyFLAP(dta,varagin)
%CLASSIFYFLAP uses activity data to classify migratory flapping flight.
%   ...

if (nargin>1)
    period = varagin{1};
else
    period = 1
end

if (nargin>2)
    plotit = varagin{2};
else
    plotit = 1
end

if (nargin>3)
    method = varagin{3};
else
    method = "kmeans"
end


  if (method == "kmeans") 
    [km,C] = kmeans(dta.act, 2); % two clusters
    threshold = mean(C);
    % dta.clust = km.cluster
  elseif (method == "hmm")
%     hmm = depmix(act ~ 1, family = poisson(), nstates = 2, 
%       data = dta[dta.act > 0, ])
%     hmmfit = fit(hmm, verbose = FALSE)
%     dta.clust = NA
%     dta.clust[dta.act > 0] = posterior(hmmfit).state
   end
  
  if (toPLOT == TRUE)
    figure; hold on;
    set(gca, 'YScale', 'log')
    histogram(dta.act(km==1));
    histogram(dta.act(km==2));
    xline(threshold)
    xlabel('Activity'); ylabel('Histogram'); legend('Low activity','High activity','Threshold')
    title(['Classification of Activity with ' method])
  end
  
  
  start = 0
  end = 0
  Duration_table = data.frame(matrix(c("2015-01-01", "2015-01-01", 
    "2015-01-01", "2015-01-01", 0, 0), nrow = 2))
  colnames(Duration_table) = c("start", "end", "Duration (h)")
  Duration_table.start = as.POSIXct(Duration_table.start, 
    tz = tz, format = "%Y-%m-%d")
  Duration_table.end = as.POSIXct(Duration_table.end, tz = tz, 
    format = "%Y-%m-%d")
  Duration_table.`Duration (h)` = as.numeric(Duration_table.`Duration (h)`)
  high_movement = as.numeric(which(table(dta.clust) == min(table(dta.clust), 
    na.rm = TRUE)))
  low_movement = as.numeric(which(table(dta.clust) == max(table(dta.clust), 
    na.rm = TRUE)))
  dta.clust[is.na(dta.clust)] = low_movement
  x = c(high_movement, low_movement, high_movement)
  idx = which(dta.clust == x[1])
  idx = idx[sapply(idx, function(i) all(dta.clust[i:(i + (length(x) - 
    1))] == x))]
  dta.clust[idx + 1] = high_movement
  x = c(low_movement, high_movement)
  start = which(dta.clust == x[1])
  start = start[sapply(start, function(i) all(dta.clust[i:(i + 
    (length(x) - 1))] == x))]
  x = c(high_movement, low_movement)
  end = which(dta.clust == x[1])
  end = end[sapply(end, function(i) all(dta.clust[i:(i + (length(x) - 
    1))] == x))]
  if (end[1] < start[1]) 
    end = end[-1]
  if (length(end) < length(start)) 
    start = start[1:length(end)]
  if (length(end) > length(start)) 
    end = end[1:length(start)]
  index = which((end - start) >= period)
  start = start[index]
  end = end[index]
  index = unlist(sapply(1:length(start), function(i) start[i]:end[i]))
  dta.clust[index] = 3
  end = c(which(dta.clust == 3)[diff(which(dta.clust == 3)) > 
    1], which(dta.clust == 3)[length(which(dta.clust == 
    3))])
  start = c(which(dta.clust == 3)[1], (which(dta.clust == 
    3)[which(diff(which(dta.clust == 3)) > 1) + 1]))
  dur = difftime(dta.date[end], dta.date[start], tz = tz, 
    units = "hours")
  info = data.frame(dta.date[start], dta.date[end], dur)
  names(info) = c("start", "end", "Duration (h)")
  Duration_table = rbind(Duration_table, info)
  if (high_movement == 1) {
    dta.clust[dta.clust == 1] = 999
    dta.clust[dta.clust == 2] = 1
    dta.clust[dta.clust == 999] = 2
  }
  Duration_table = Duration_table[-c(1, 2), ]
  dta.clust[dta.act == 0] = 0
  return(list(timetable = Duration_table, classification = dta.clust, 
    no_movement = 0, low_movement = 1, high_movement = 2, 
    migration = 3, threshold = threshold))
}


end

