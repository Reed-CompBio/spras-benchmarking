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

export function extractDatasetType(name: string): { type: string, name: string } {
    let newType;
    let newName;
    
    for (let type of dataTypes) {
        if (name.startsWith(type)) {
            newType = type;
            newName = name.substring(type.length);
            break;
        }
    }

    // We add the !newName there for type-checking purposes.
    if (!newType || !newName) {
        throw new Error(`Dataset name should begin with a type (one of ${dataTypes})`)
    }

    return { 
        type: newType,
        name: newName
    }
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


    // We didn't get a data type in the first passthrough - lets extract the data
    // type from the name
    if (!dataType) {
        const { type, name } = extractDatasetType(datasetName);
        dataType = type;
        datasetName = name;
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
