import os

# Настройки
OUTPUT_FILE = "project_dump.txt"
# Расширения файлов, которые мы хотим включить в дамп
INCLUDE_EXTENSIONS = {'.py', '.json', '.yml', '.txt', '.md'}
# Папки, которые нужно игнорировать (чтобы не забить файл мусором)
IGNORE_DIRS = {'.git', '__pycache__', 'data', 'optimization_results', 'logs', '.venv', 'venv'}

def export_project_code():
    project_root = os.getcwd()
    
    with open(OUTPUT_FILE, "w", encoding="utf-8") as outfile:
        outfile.write(f"PROJECT DUMP: {project_root}\n")
        outfile.write("="*60 + "\n\n")
        
        for root, dirs, files in os.walk(project_root):
            # Модифицируем dirs на месте, чтобы os.walk не заходил в игнорируемые папки
            dirs[:] = [d for d in dirs if d not in IGNORE_DIRS]
            
            for file in files:
                # Игнорируем сам файл дампа и скрытые файлы (кроме .json/.yml/.md)
                if file == OUTPUT_FILE or file.startswith('.'):
                    continue
                
                file_path = os.path.join(root, file)
                extension = os.path.splitext(file)[1]
                
                if extension in INCLUDE_EXTENSIONS:
                    relative_path = os.path.relpath(file_path, project_root)
                    
                    try:
                        with open(file_path, "r", encoding="utf-8") as infile:
                            content = infile.read()
                            
                        outfile.write(f"\n{'='*20} FILE: {relative_path} {'='*20}\n")
                        outfile.write(content)
                        outfile.write(f"\n{'='*20} END OF FILE: {relative_path} {'='*20}\n")
                        print(f"Added: {relative_path}")
                        
                    except Exception as e:
                        print(f"Skipped {relative_path}: {e}")

    print(f"\nDone! All code exported to {OUTPUT_FILE}")

if __name__ == "__main__":
    export_project_code()
