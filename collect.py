#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from pathlib import Path

# ==============================================================================
# CONFIGURACIÓN DE CARPETAS Y FILTROS
# ==============================================================================
# Define aquí los nombres de las carpetas específicas que deseas consolidar.
# Por ejemplo: TARGET_FOLDERS = ['src', 'modulos', 'test']
# Si dejas la lista vacía [], el script escaneará todo el proyecto desde la raíz.
TARGET_FOLDERS = ['src', 'tests', 'Theory', 'pages']
TARGET_FILES = ['CMakeLists.txt', 'my_project.md']

# Nombre del archivo unificado resultante
OUTPUT_FILENAME = "all_fortran.txt"

def clean_path(path, root_dir):
    """Devuelve la ruta relativa para una visualización limpia en la IA."""
    try:
        return path.relative_to(root_dir)
    except ValueError:
        return path


def collect_files(root_dir):
    """Recorre únicamente las carpetas seleccionadas por el usuario."""
    valid_files = []
    
    # Determinar las rutas a escanear
    folders_to_scan = []
    if TARGET_FOLDERS:
        for folder in TARGET_FOLDERS:
            folder_path = root_dir / folder
            if folder_path.exists() and folder_path.is_dir():
                folders_to_scan.append(folder_path)
            else:
                print(f"⚠️ Advertencia: La carpeta configurada '{folder}' no existe en este directorio.")
    else:
        # Si no hay carpetas definidas, se escanea el directorio raíz completo
        folders_to_scan = [root_dir]

    # Escaneo de archivos en las rutas elegidas
    for scan_path in folders_to_scan:
        for dirpath, _, filenames in os.walk(scan_path):
            for filename in filenames:
                file_path = Path(dirpath) / filename
                valid_files.append(file_path)

    if 'TARGET_FILES' in globals():
        for filename in TARGET_FILES:
            file_path = root_dir / filename
            if file_path.exists() and file_path.is_file():
                if file_path not in valid_files: # Evitar duplicados
                    valid_files.append(file_path)
                                
    # Ordenar alfabéticamente para mantener una estructura predecible
    valid_files.sort()
    return valid_files

def build_consolidated_file(root_path, files):
    """Genera el gran archivo de texto plano estructurando el contenido."""
    output_path = root_path / OUTPUT_FILENAME
    
    print(f"📝 Creando archivo unificado en: {output_path.name}")
    
    try:
        with open(output_path, 'w', encoding='utf-8', errors='replace') as outfile:
            # Cabecera estructurada
            outfile.write("===================================================\n")
            outfile.write("CONSOLIDADO DE CÓDIGO FUENTE DE PROYECTO FORTRAN\n")
            outfile.write("===================================================\n\n")
            outfile.write("Este archivo contiene la estructura completa del código del proyecto\n")
            
            # Sección 1: Índice
            outfile.write("--- 1. ÍNDICE DE ARCHIVOS INCLUIDOS ---\n")
            for i, file_path in enumerate(files, 1):
                rel_path = clean_path(file_path, root_path)
                outfile.write(f"[{i}] {rel_path}\n")
            outfile.write("\n" + "="*50 + "\n\n")
            
            # Sección 2: Contenido detallado de cada archivo
            outfile.write("--- 2. CONTENIDO DE LOS ARCHIVOS ---\n\n")
            
            for file_path in files:
                rel_path = clean_path(file_path, root_path)
                print(f" -> Procesando: {rel_path}")
                
                outfile.write(f"#### INICIO ARCHIVO: {rel_path} ####\n")
                outfile.write(f"Ruta: {rel_path}\n")
                outfile.write("-" * 40 + "\n")
                
                try:
                    with open(file_path, 'r', encoding='utf-8', errors='replace') as infile:
                        content = infile.read()
                        outfile.write(content)
                except Exception as e:
                    outfile.write(f"[ERROR al leer este archivo: {str(e)}]\n")
                
                outfile.write("\n")
                outfile.write(f"#### FIN ARCHIVO: {rel_path} ####\n")
                outfile.write("=" * 50 + "\n\n")
                
        print(f"✨ ¡Éxito! Archivo '{OUTPUT_FILENAME}' generado con {len(files)} archivos unificados.")
        
    except Exception as e:
        print(f"❌ Error al escribir el archivo unificado: {e}", file=sys.stderr)

def main():
    current_directory = Path('.').resolve()
    
    print("==============================================")
    print("   Consolidador de Código Fortran             ")
    print("==============================================")
    print(f"Carpeta base del proyecto: {current_directory}")
    if TARGET_FOLDERS:
        print(f"Carpetas seleccionadas a escanear: {', '.join(TARGET_FOLDERS)}")
    else:
        print("Escaneando: Todo el proyecto (raíz y subcarpetas)")
    
    # Recolectar archivos válidos
    files_to_process = collect_files(current_directory)
    
    if not files_to_process:
        print("⚠️ No se encontraron archivos válidos en las carpetas seleccionadas.")
        return
        
    print(f"🔍 Encontrados {len(files_to_process)} archivos aptos para consolidación.")
    
    # Construir el archivo final unificado
    build_consolidated_file(current_directory, files_to_process)

if __name__ == "__main__":
    main()