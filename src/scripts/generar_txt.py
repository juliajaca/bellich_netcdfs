
from netCDF4 import Dataset

def generar_txt(nombre_nc, txt_filename):

    # nombre_nc = nombre_nc.replace(":", "_")
    ncfile = Dataset(nombre_nc, "r")
    # txt_filename = "albujon_flow_hourly_output.txt"

    with open(txt_filename, "w") as f:
        # Escribir el formato
        f.write(f"Format:\n\t{ncfile.file_format}\n")

        # Escribir atributos globales
        f.write("\nGlobal Attributes:\n")
        for attr in ncfile.ncattrs():
            f.write(f"\t{attr:<15} = '{getattr(ncfile, attr)}'\n")

        # Escribir dimensiones
        f.write("\nDimensions:\n")
        for dim_name, dim in ncfile.dimensions.items():
            f.write(f"\t{dim_name:<15} = {len(dim)}\n")

        # Escribir variables
        f.write("\nVariables:\n")
        for var_name, var in ncfile.variables.items():
            # Obtener dimensiones y tipo de datos
            dimensions = ", ".join(var.dimensions)
            dtype = var.dtype

            # Escribir información de la variable
            f.write(f"    {var_name:<15}\n")
            f.write(f"\tSize: {var.shape}\n")
            f.write(f"\tDimensions: {dimensions}\n")
            f.write(f"\tDatatype: {dtype}\n")
            
            # Escribir atributos de la variable
            if var.ncattrs():
                f.write("\tAttributes:\n")
                for attr in var.ncattrs():
                    f.write(f"\t\t{attr:<15} = '{getattr(var, attr)}'\n")

    # Cerrar el archivo NetCDF
    ncfile.close()

    print(f"Archivo '{txt_filename}' generado con éxito.")
