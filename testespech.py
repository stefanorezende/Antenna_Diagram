import speech_recognition as sr

def ouvir_microfone():
    # Cria um objeto Recognizer
    recognizer = sr.Recognizer()

    # Abre o microfone para captura
    with sr.Microphone() as source:
        print("Diga alguma coisa...")
        recognizer.adjust_for_ambient_noise(source)  # Ajusta para o ruído ambiente
        audio = recognizer.listen(source)  # Captura o áudio do microfone

    try:
        # Usa o recognizer para converter áudio em texto
        texto = recognizer.recognize_google(audio, language='pt-BR')
        print("Você disse:", texto)
    except sr.UnknownValueError:
        print("Não foi possível entender o áudio")
    except sr.RequestError as e:
        print("Erro ao chamar o serviço de reconhecimento de fala do Google; {0}".format(e))

# Chama a função para começar a ouvir
ouvir_microfone()
