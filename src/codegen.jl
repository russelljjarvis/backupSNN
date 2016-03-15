export @codegen
macro codegen()
  str="""
	function print(ex...)
		println(ex...)
	end"""
  eval(parse(str))
end
