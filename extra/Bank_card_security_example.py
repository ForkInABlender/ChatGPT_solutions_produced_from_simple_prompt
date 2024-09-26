# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""
This is what is written on cards and how they're verified.

"""

from Crypto.Cipher import AES
from Crypto.PublicKey import RSA
from Crypto.Cipher import PKCS1_OAEP
from Crypto.Random import get_random_bytes
from Crypto.Hash import SHA256, HMAC
from Crypto.Signature import pkcs1_15
import base64

# Generate RSA keys (This would typically be done inside an HSM)
def generate_rsa_keys():
    key = RSA.generate(2048)
    private_key = key.export_key()
    public_key = key.publickey().export_key()
    return private_key, public_key

# AES encryption for credit card data
def aes_encrypt(key, data):
    cipher = AES.new(key, AES.MODE_GCM)
    ciphertext, tag = cipher.encrypt_and_digest(data)
    return base64.b64encode(cipher.nonce + tag + ciphertext)

# AES decryption for credit card data
def aes_decrypt(key, encrypted_data):
    decoded = base64.b64decode(encrypted_data)
    nonce = decoded[:16]
    tag = decoded[16:32]
    ciphertext = decoded[32:]
    cipher = AES.new(key, AES.MODE_GCM, nonce=nonce)
    return cipher.decrypt_and_verify(ciphertext, tag)

# RSA for Key Management (to securely store the AES key)
def rsa_encrypt_key(public_key, key):
    rsa_public_key = RSA.import_key(public_key)
    cipher = PKCS1_OAEP.new(rsa_public_key)
    encrypted_key = cipher.encrypt(key)
    return base64.b64encode(encrypted_key)

def rsa_decrypt_key(private_key, encrypted_key):
    rsa_private_key = RSA.import_key(private_key)
    cipher = PKCS1_OAEP.new(rsa_private_key)
    decrypted_key = cipher.decrypt(base64.b64decode(encrypted_key))
    return decrypted_key

# HMAC for Integrity Check
def generate_hmac(key, data):
    hmac_instance = HMAC.new(key, digestmod=SHA256)
    hmac_instance.update(data)
    return hmac_instance.hexdigest()

def verify_hmac(key, data, hmac_to_verify):
    hmac_instance = HMAC.new(key, digestmod=SHA256)
    hmac_instance.update(data)
    try:
        hmac_instance.hexverify(hmac_to_verify)
        return True
    except ValueError:
        return False

# Example Usage for Credit Card Security

# Example credit card details (PAN + other sensitive data)
card_data = b'4111111111111111|12/24|123'  # Example PAN + Expiry + CVV

# Generate a 256-bit AES key
aes_key = get_random_bytes(32)

# Encrypt the card data
encrypted_card_data = aes_encrypt(aes_key, card_data)
print("Encrypted Card Data:", encrypted_card_data)

# Decrypt the card data
decrypted_card_data = aes_decrypt(aes_key, encrypted_card_data)
print("Decrypted Card Data:", decrypted_card_data)

# RSA Key Generation (for managing AES keys securely)
private_key, public_key = generate_rsa_keys()

# Encrypt the AES key using RSA (to mimic key management in HSM)
encrypted_aes_key = rsa_encrypt_key(public_key, aes_key)
print("Encrypted AES Key:", encrypted_aes_key)

# Decrypt the AES key using RSA
decrypted_aes_key = rsa_decrypt_key(private_key, encrypted_aes_key)
print("Decrypted AES Key:", decrypted_aes_key)

# HMAC for ensuring integrity of the card data
hmac_value = generate_hmac(aes_key, card_data)
print("Generated HMAC:", hmac_value)

# Verify HMAC
is_valid_hmac = verify_hmac(aes_key, card_data, hmac_value)
print("HMAC Verified:", is_valid_hmac)
