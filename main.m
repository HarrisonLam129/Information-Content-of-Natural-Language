a = readmatrix('II-19-2-dataA.txt');
a = a(1:400,:);

prob = zeros(1,27);
logprob = zeros(1,27);
for i = 1:27
    prob(i) = sum(a(:)==i-1)/numel(a);
    if prob(i) > 0
        logprob(i) = log(prob(i))/log(2);
    else
        %If letter does not appear in the sample, set logp=0 so plogp=0
        logprob(i) = 0;
    end
end
entropy = -dot(prob, logprob);

%Huffman
huffmancws = strings(1, 27);
tempindex = 1:27;
probtemp = prob;
for i = 1:26
    %Find the 2 groups with the lowest probabilities
    [minprob1, index1] = min(probtemp);
    probtemp(index1) = probtemp(index1)+2;
    [minprob2, index2] = min(probtemp);
    append1 = tempindex(1,index1,:);
    append2 = tempindex(1,index2,:);
    %Add 0s to all the elements in the lowest group
    for j = 1:length(append2)
        if append2(j) ~= 0
            huffmancws(append2(j)) = append('0', huffmancws(append2(j)));
        end
    end
    %Add 1s to all the elements in the 2nd lowest group
    for j = 1:length(append1)
        if append1(j) ~= 0
            huffmancws(append1(j)) = append('1', huffmancws(append1(j)));
        end
    end
    %Combine the groups
    probtemp(index2) = probtemp(index2) + probtemp(index1)-2;
    probtemp(index1) = 2;
    for k = 1:length(append1)
        if append1(k) ~= 0
            if isempty(find(append2 == 0, 1))
                tempindex(1,index2,length(append2)+k) = append1(k);
            else
                tempindex(1,index2,find(append2 == 0, 1)-1+k) = append1(k);
            end
        end
    end
    tempindex(1,index1,:) = 0;
end
el1 = dot(strlength(huffmancws),prob);

%Shannon-Fano
sfcws = strings(1, 27);
sflength = ceil(-logprob);
%Order the lengths from lowest to highest
[sflength, orderstat] = sort(sflength);
for i = 1:27
    %Loop through all the words of desired length
    for k = 1:2^sflength(i)
        flag = 0;
        bin = dec2bin(k-1,sflength(i));
        %Check if the current word is such that the code is prefix-free
        for j = 1:i-1
            str = sfcws(orderstat(j));
            if bin(1:strlength(str)) == str
                flag = 1;
                break
            end
        end
        if flag == 0
            sfcws(orderstat(i)) = bin;
            break
        end
    end
end
el2 = dot(strlength(sfcws),prob);

b = a';
b = b(:);


%Counting #times each letter appears right after letterIndex appears
letterIndex = 9;
afterletter = [];
for i = 1:9999
    if b(i) == letterIndex
        afterletter = [afterletter b(i+1)];
    end
end
prob2 = zeros(1,27);
for i = 1:27
    prob2(i) = sum(afterletter(:)==i-1)/numel(afterletter);
end
plot(0:26,prob);
hold on
plot(0:26,prob2)
ylabel('Probability')
xlabel('Letter');
legend('General', append('After letter ', num2str(letterIndex)));

%HuffmanPair
pairprob = zeros(1,729);
a = reshape(a',2,5000);
for i = 1:27
    for j = 1:27
        paircount = (a==[i-1;j-1]);
        %Using AND operator to check if both characters match
        paircount = sum(paircount(1,:)&paircount(2,:));
        if paircount > 0
            pairprob(27*(i-1)+j) = paircount/5000;
        else
            %Estimate 0 probability with Bernoulli probability
            pairprob(27*(i-1)+j) = prob(i)*prob(j);
        end
    end
end
%Rescaling probabilities
pairprob = pairprob/sum(pairprob);

huffmanpaircws = strings(1, 729);
tempindex = 1:729;
probtemp = pairprob;
for i = 1:728
    %Find the 2 groups with the lowest probabilities
    [minprob1, index1] = min(probtemp);
    probtemp(index1) = probtemp(index1)+2;
    [minprob2, index2] = min(probtemp);
    append1 = tempindex(1,index1,:);
    append2 = tempindex(1,index2,:);
    %Add 0s to all the elements in the lowest group
    for j = 1:length(append2)
        if append2(j) ~= 0
            huffmanpaircws(append2(j)) = append('0', huffmanpaircws(append2(j)));
        end
    end
    %Add 1s to all the elements in the 2nd lowest group
    for j = 1:length(append1)
        if append1(j) ~= 0
            huffmanpaircws(append1(j)) = append('1', huffmanpaircws(append1(j)));
        end
    end
    %Combine the groups
    probtemp(index2) = probtemp(index2) + probtemp(index1)-2;
    probtemp(index1) = inf;
    for k = 1:length(append1)
        if append1(k) ~= 0
            if isempty(find(append2 == 0, 1))
                tempindex(1,index2,length(append2)+k) = append1(k);
            else
                tempindex(1,index2,find(append2 == 0, 1)-1+k) = append1(k);
            end
        end
    end
    tempindex(1,index1,:) = 0;
end
el3 = dot(strlength(huffmanpaircws), pairprob);



%Random generated text
cprob = cumsum([0, prob]);
[~, ~, randomtext] = histcounts(rand(1,10000),cprob);
randomtext = randomtext-1;
randomtext = reshape(randomtext,2,5000);
randomtextcount = zeros(1,27);
randomtextpaircount = zeros(1,729);
bernoullipairprob = zeros(1,729);
for i = 1:27
    randomtextcount(i) = sum(randomtext(:)==i-1);
    for j = 1:27
        paircount = (randomtext==[i-1;j-1]);
        randomtextpaircount(27*(i-1)+j) = sum(paircount(1,:)&paircount(2,:));
        bernoullipairprob(27*(i-1)+j) = prob(i)*prob(j);
    end
end

%HuffmanPairBernoulli
huffmanpaircws2 = strings(1, 729);
tempindex = 1:729;
probtemp = bernoullipairprob;
for i = 1:728
    %Find the 2 groups with the lowest probabilities
    [minprob1, index1] = min(probtemp);
    probtemp(index1) = probtemp(index1)+2;
    [minprob2, index2] = min(probtemp);
    append1 = tempindex(1,index1,:);
    append2 = tempindex(1,index2,:);
    %Add 0s to all the elements in the lowest group
    for j = 1:length(append2)
        if append2(j) ~= 0
            huffmanpaircws2(append2(j)) = append('0', huffmanpaircws2(append2(j)));
        end
    end
    %Add 1s to all the elements in the 2nd lowest group
    for j = 1:length(append1)
        if append1(j) ~= 0
            huffmanpaircws2(append1(j)) = append('1', huffmanpaircws2(append1(j)));
        end
    end
    %Combine the groups
    probtemp(index2) = probtemp(index2) + probtemp(index1)-2;
    probtemp(index1) = inf;
    for k = 1:length(append1)
        if append1(k) ~= 0
            if isempty(find(append2 == 0, 1))
                tempindex(1,index2,length(append2)+k) = append1(k);
            else
                tempindex(1,index2,find(append2 == 0, 1)-1+k) = append1(k);
            end
        end
    end
    tempindex(1,index1,:) = 0;
end
singlelength = dot(randomtextcount, strlength(huffmancws))
pairlength = dot(randomtextpaircount, strlength(huffmanpaircws2))