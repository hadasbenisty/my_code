function [labs, numlabs] = getModuleLabels(module,  mainLabels)
labs = cell(length(module), 1);
for k=1:length(module)
    for level = length(module{k}):-1:1
        if any(strcmp(module{k}{level}, mainLabels))
            labs{k} = module{k}{level};
numlabs(k) = find(strcmp(module{k}{level}, mainLabels));
            break;
        end
    end
end
