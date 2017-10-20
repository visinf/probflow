function [val] = getParam(paramStruct,param,defaultVal)

    if isfield(paramStruct,param)
            val = paramStruct.(param);
        else
            val = defaultVal;
    end
end

