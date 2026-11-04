import java.nio.charset.StandardCharsets
import java.nio.file.Files

/**
 * Expands environment variables in CSV inputs used by NFTest configurations.
 */
class NFTestCsvTools {
    /**
     * Expand $NAME and ${NAME} references using the environment inherited by
     * Nextflow. Undefined variables fail instead of leaving unusable paths.
     */
    String expandEnvironmentVariables(String content) {
        if (content == null) {
            return null
        }

        def variablePattern = ~/\$(?:\{([A-Za-z_][A-Za-z0-9_]*)\}|([A-Za-z_][A-Za-z0-9_]*))/

        return content.replaceAll(variablePattern) { match, bracedName, bareName ->
            def variableName = bracedName ?: bareName
            def variableValue = System.getenv(variableName)

            if (variableValue == null) {
                throw new IllegalArgumentException(
                    "Environment variable ${variableName} referenced by test input is not defined"
                )
            }

            return variableValue
        }
    }

    /** Check that every input path referenced by a test CSV exists. */
    void requireInputFiles(Object csvPath, String pathColumn, boolean allowDirectory = false) {
        def csvFile = new File(csvPath.toString())
        def lines = Files.readAllLines(csvFile.toPath(), StandardCharsets.UTF_8)
        if (!lines) {
            throw new IllegalArgumentException("Empty test input CSV: ${csvFile}")
        }
        def columns = lines.first().split(',', -1)*.trim()
        def columnIndex = columns.indexOf(pathColumn)
        if (columnIndex < 0) {
            throw new IllegalArgumentException("${csvFile} is missing column ${pathColumn}")
        }
        lines.drop(1).eachWithIndex { line, rowIndex ->
            if (!line.trim()) return
            def fields = line.split(',', -1)
            def value = columnIndex < fields.size() ? fields[columnIndex].trim() : ''
            if (!value) {
                throw new IllegalArgumentException(
                    "Empty ${pathColumn} in ${csvFile} at row ${rowIndex + 2}")
            }
            def input = new File(value)
            if (!input.isAbsolute()) input = new File(csvFile.parentFile, value)
            if (!(input.isFile() || (allowDirectory && input.isDirectory()))) {
                throw new IllegalArgumentException(
                    "Missing ${pathColumn} input from ${csvFile} row ${rowIndex + 2}: ${input}")
            }
        }
    }

    /**
     * Expand environment variables in a CSV and write the result beneath the
     * Nextflow work directory. NFTest removes this with the test work tree.
     */
    String prepareInputCsv(Object sourceCsv, String testName, Object launchDir) {
        def sourceFile = new File(sourceCsv.toString())
        if (!sourceFile.isFile()) {
            throw new IllegalArgumentException("Test input CSV does not exist: ${sourceFile}")
        }

        def workRoot = new File(System.getenv('NXF_WORK') ?: "${launchDir}/work")
        def inputDirectory = new File(workRoot, 'nftest-input')
        Files.createDirectories(inputDirectory.toPath())

        def safeTestName = testName.replaceAll(/[^A-Za-z0-9_.-]/, '-')
        def expandedCsv = File.createTempFile("${safeTestName}-", '.csv', inputDirectory)
        def sourceContent = Files.readString(sourceFile.toPath(), StandardCharsets.UTF_8)
        def expandedContent = expandEnvironmentVariables(sourceContent)

        Files.writeString(expandedCsv.toPath(), expandedContent, StandardCharsets.UTF_8)
        expandedCsv.deleteOnExit()

        return expandedCsv.absolutePath
    }
}

return new NFTestCsvTools()
