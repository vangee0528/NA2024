# 本程序用于一键生成文件结构树
# 在指定文件夹下运行本程序，即可生成文件结构树,并将结果复制到剪贴板。
import os
import pyperclip

# 在此处填写需要忽略的文件扩展名
ignore_extensions = {'.json', '.exe' ,'.gitkeep'}

# 添加一个开关来控制是否忽略以`.`开头的文件夹
ignore_dot_folders = True  

def generate_file_structure_tree(folder_path, prefix=''):
    files = []
    dirs = []
    
    items = os.listdir(folder_path)
    
    for item in items:
        path = os.path.join(folder_path, item)
        
        if os.path.isfile(path) and os.path.splitext(item)[1] in ignore_extensions:
            continue
        
        if os.path.isdir(path):
            if ignore_dot_folders and item.startswith('.'):
                continue
            dirs.append(item)
        else:
            files.append(item)

    files.sort()
    dirs.sort()
    
    sorted_items = files + dirs

    tree = []
    for i, item in enumerate(sorted_items):
        path = os.path.join(folder_path, item)
        is_last = i == len(sorted_items) - 1
        
        if os.path.isdir(path):
            tree.append(f"{prefix}{'└── ' if is_last else '├── '}{item}/")
            tree.extend(generate_file_structure_tree(path, prefix + ('    ' if is_last else '│   ')))
        else:
            tree.append(f"{prefix}{'└── ' if is_last else '├── '}{item}")

    return tree

def format_tree(tree):
    current_folder_name = os.path.basename(os.getcwd())  # 获取当前文件夹名称
    formatted_tree = [f"{current_folder_name}/"]
    for line in tree:
        formatted_tree.append("    " + line)
    return "\n".join(formatted_tree)

def main():
    folder_path = './' # 默认为当前文件夹

    tree = generate_file_structure_tree(folder_path)
    formatted_output = format_tree(tree)
        
    pyperclip.copy(formatted_output)
        
    print("\n文件结构：")
    print(formatted_output)
    print("\n已复制到剪贴板。")

if __name__ == "__main__":
    main()
