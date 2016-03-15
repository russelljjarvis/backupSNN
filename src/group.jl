export makeGroup,connect,@group
function makeGroup(groupName,groupSize)
  group=Dict()
  group["N"]=sum(groupSize)
  startIndex=cumsum(groupSize)-groupSize+1
  for (idx,name) in enumerate(groupName)
    group[name]=startIndex[idx]:startIndex[idx]+groupSize[idx]-1
  end
  return group
end

macro connect(groups,preGroupName,postGroupName,condition,cmd)
  var=cmd.args[1].args[1]
  lhs=cmd.args[1]
  rhs=cmd.args[2]
  if !isdefined(var)
    eval(:($var=zeros($groups["N"],$groups["N"])))
  end
  ex=quote
    for post=1:length($groups[$postGroupName])
      for pre=1:length($groups[$preGroupName])
        if $condition
          rhsValue=$rhs
          pre=pre+$(groups)[$preGroupName][1]-1
          post=post+$(groups)[$postGroupName][1]-1
          $lhs=rhsValue
        end
      end
    end
  end
  eval(ex)
end

macro group(groups,groupName,cmd)
  var=cmd.args[1].args[1]
  rhs=cmd.args[2]
  lhs=cmd.args[1]
  if !isdefined(var)
    eval(:($var=zeros($groups["N"])))
  end
  ex=quote
    for i=1:length($groups[$groupName])
      rhsValue=$rhs
      i=i+$groups[$groupName][1]-1
      $lhs=rhsValue
    end
  end
  eval(ex)
end
