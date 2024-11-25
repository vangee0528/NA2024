def decimal_to_binary(num, precision=10):
    # 提取符号位
    sign = '-' if num < 0 else ''
    num = abs(num)
    
    # 整数部分转换
    integer_part = int(num)
    fractional_part = num - integer_part
    
    binary_integer = bin(integer_part)[2:]  # 整数部分直接使用 bin() 转换
    
    # 小数部分转换
    binary_fraction = []
    for _ in range(precision):
        fractional_part *= 2
        bit = int(fractional_part)
        binary_fraction.append(str(bit))
        fractional_part -= bit
        if fractional_part == 0:  # 精确转换完成
            break
    
    binary_fraction_str = ''.join(binary_fraction)
    
    return f"{sign}{binary_integer}.{binary_fraction_str}"


# 测试用例
number = 477
precision = 20  # 保留20位小数部分
binary_result = decimal_to_binary(number, precision)
print(f"({number})_10 = ({binary_result})_2")
