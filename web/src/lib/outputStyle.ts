interface Output {
    datasetName: string;
    algorithm: string;
    paramsHash: string;
}

export function parseOutputString(str: string): Output {
    const components = str.split("-");
    let datasetName;
    let algorithm;
    let paramsHash;

    if (components.length === 4) {
        // We ignore -params-
        [datasetName, algorithm, , paramsHash] = components
    } else if (components.length === 3) {
        [datasetName, algorithm, paramsHash] = components
    } else {
        throw new Error(`Unexpected length of components in ${components}.`)
    }

    return {
        datasetName,
        algorithm,
        paramsHash
    }
}

export function asFolderName(output: Output): string {
    return `${output.datasetName}-${output.algorithm}-params-${output.paramsHash}`
}
