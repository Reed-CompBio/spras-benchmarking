const dataTypes = [
    'pra',
    'dmmm'
]

interface Output {
    dataType: string;
    datasetName: string;
    algorithm: string;
    paramsHash: string;
}

export function parseOutputString(str: string): Output {
    const components = str.split("-");
    let dataType;
    let datasetName;
    let algorithm;
    let paramsHash;

    if (components.length === 4) {
        if (dataTypes.includes(components[0])) {
            // This is a slug URL (type-...)
            [dataType, datasetName, algorithm, paramsHash] = components
        } else {
            // This is fetched straight from the folder - we ignore -params-
            [datasetName, algorithm, , paramsHash] = components
        }
    } else if (components.length === 3) {
        [datasetName, algorithm, paramsHash] = components
    } else {
        throw new Error(`Unexpected length of components in ${components}.`)
    }


    // We didn't get a data type in the first passthrough - lets see if the dataset name
    // starts with a type.
    if (!dataType) {
        for (let type of dataTypes) {
            if (datasetName.startsWith(type)) {
                dataType = type;
                datasetName = datasetName.substring(type.length);
                break;
            }
        }
    }

    if (!dataType) {
        throw new Error(`Dataset name should begin with a type (one of ${dataTypes})`)
    }

    return {
        dataType,
        datasetName,
        algorithm,
        paramsHash
    }
}

export function styleOutput(output: Output): string {
    return `${output.dataType}-${output.datasetName}-${output.algorithm}-${output.paramsHash}`
}

export function asFolderName(output: Output): string {
    return `${output.dataType}${output.datasetName}-${output.algorithm}-params-${output.paramsHash}`
}
