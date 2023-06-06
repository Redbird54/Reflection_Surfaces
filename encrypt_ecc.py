from tinyec.ec import SubGroup, Curve
from Crypto.Cipher import AES
import hashlib, secrets, binascii, struct
import sys

def ecc_point_to_256_bit_key(point):
    sha = hashlib.sha256(int.to_bytes(point.x, 32, 'big'))
    sha.update(int.to_bytes(point.y, 32, 'big'))
    return sha.digest()

p = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f
n = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
a = 0x0000000000000000000000000000000000000000000000000000000000000000
b = 0x0000000000000000000000000000000000000000000000000000000000000007
g = (0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
    0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8)
h = 1
curve = Curve(a, b, SubGroup(p, g, n, h), 'secp256k1')
privKey = secrets.randbelow(curve.field.n)
pubKey = privKey * curve.g

def encrypt (msg1, name):
    msg = msg1.encode('UTF-8')

    ciphertextPrivKey = secrets.randbelow(curve.field.n)
    secretKey = ecc_point_to_256_bit_key(ciphertextPrivKey * pubKey)
    ciphertextPubKey = ciphertextPrivKey * curve.g
    aesCipherEnc = AES.new(secretKey, AES.MODE_GCM)
    ciphertext, authTag = aesCipherEnc.encrypt_and_digest(msg)
    encryptedMsg = (ciphertext, aesCipherEnc.nonce, authTag, ciphertextPubKey)

    print(name)
    print("Encrypted msg:", float(int(str(binascii.hexlify(encryptedMsg[0]))[2:-1], 16)))
    return encryptedMsg

def decrypt (encryptedMsg):
    (ciphertext, nonce, authTag, ciphertextPubKey) = encryptedMsg
    secretKey = ecc_point_to_256_bit_key(privKey * ciphertextPubKey)
    aesCipherDec = AES.new(secretKey, AES.MODE_GCM, nonce)
    decryptedMsg = aesCipherDec.decrypt_and_verify(ciphertext, authTag)


    print("Decrypted msg:", str(decryptedMsg)[2:-1])
    return float(decryptedMsg)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise Exception("Must have 2 arguments")
        sys.exit()
    msg = encrypt(sys.argv[1], str(sys.argv[2]))
    output = decrypt(msg)