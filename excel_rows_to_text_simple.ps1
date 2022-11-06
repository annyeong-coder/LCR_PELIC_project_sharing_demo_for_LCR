$counter = 1

Get-Content "answer_ID_text_simplified_v05full.txt" | foreach {
    Set-Content -Path "$counter.txt" -Value $_
    $counter++
}