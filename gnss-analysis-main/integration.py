if __name__ == '__main__':
    try:
        if len(sys.argv) < 2:
            print("Usage: python script_name.py <file_path>")
            sys.exit(1)

        file_path = sys.argv[1]
        print(f"Processing file: {file_path}")  # Debug: output the file path being processed

        ecef_List = process(file_path)
        csv_FromList = pd.DataFrame(ecef_List)
        csv_FromList.to_csv("Measurements_Output.csv", index=False)

        file_measurements_output = 'android_position.csv'
        csv_to_kml(file_measurements_output, 'output_kml.kml')

    except Exception as e:
        print(f"An error occurred: {e}")
