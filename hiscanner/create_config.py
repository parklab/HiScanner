import json
import os

def _check_file(fn):
    '''
    Checks if the file exists and is not empty. If not, return True.
    '''
    if not os.path.exists(fn):
        return True
    if os.stat(fn).st_size == 0:
        return True
    return False
def main():
    print("You will be prompted to enter file paths and parameters required for the analysis.")
    print("If you prefer to use the default values, simply press Enter when prompted.")
    print("-------------------------------------------------------------")
    proceed = input("Do you wish to proceed with creating a JSON configuration file? (y/n): ").strip().lower()
    if proceed != "y":
        print("Initialization cancelled. Exiting program.")
        exit()

    print("\nPlease enter the required information below:")

    data = {}
    fields = {
        "gatk_vcf": "Enter the path of the VCF containing raw GATK SNP calls: ",
        "phase_file": "Enter the path of the VCF containing phased germline SNPs: ",
        "bin_path": "Enter the path to the folder containing BIC-seq2 norm output: "
    }

    for field, prompt in fields.items():
        while True:
            data[field] = input(prompt)
            if not _check_file(data[field]):
                break
            print("File does not exist or is empty. Please re-enter.")

    json_msg = "Enter the path where the output JSON file should be saved: "
    wgd_msg = "Enter the MAX_WGD value (integer, default 1): "
    germ_msg = "Enter the name of the bulk sample: "
    sc_msg = "Enter the single cell names (if more than one cell, separate by commas: "
    stem_msg = "Enter the path for the output folder: "
    lambda_msg = "Enter the LAMBDA value (integer, default 100): "
    j_msg = "Enter the number of jobs to run in parallel (integer, default 8): "


    data["germline"] = input(germ_msg)
    data["singlecell"] = input(sc_msg)
    data["stem"] = input(stem_msg)

    # create the stem directory if it doesn't exist
    os.makedirs(data["stem"], exist_ok=True)

    data["MAX_WGD"] = input(wgd_msg)
    if data["MAX_WGD"] == '':
        data["MAX_WGD"] = 1
    data['LAMBDA'] = input(lambda_msg)
    if data['LAMBDA'] == '':
        data['LAMBDA'] = 100
    data['j'] = input(j_msg)
    if data['j'] == '':
        data['j'] = 8
    
    output_file = input(json_msg)
    while not output_file.endswith(".json"):
        print("Output file must be a JSON file. Please re-enter.")
        output_file = input(json_msg)

    while int(data["j"]) <= 0:
        print("Number of jobs must be a positive integer. Please re-enter.")
        data['j'] = input(j_msg)

    while int(data["MAX_WGD"]) <= 0:
        print("MAX_WGD must be a positive integer. Please re-enter.")
        data["MAX_WGD"] = input(wgd_msg)

    while int(data["LAMBDA"]) <= 0:
        print("LAMBDA must be a positive integer. Please re-enter.")
        data['LAMBDA'] = input(lambda_msg)

    while not os.path.exists(data["stem"]):
        print("Output folder does not exist or cannot be created. Please re-enter.")
        data["stem"] = input(stem_msg)

    with open(output_file, 'w') as outfile:
        json.dump(data, outfile, indent=4)

    print(f"Configuration data has been successfully written to {output_file}")

if __name__ == "__main__":
    main()
